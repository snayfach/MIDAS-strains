#!/usr/bin/env python

import os, csv, gzip, sys, itertools, numpy as np
from collections import defaultdict
from scipy.optimize import nnls

def parse_args():
	import argparse
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""
description: use reads mapped to pangenomes to estimate strain frequencies

usage: midas_strains.py --indir <indir> [options]""",
		epilog="""Examples: 1) ...""")
	parser.add_argument('--indir', metavar='CHAR', required=True,
		help="""Path to input directory with MIDAS output for one sample
Should contain subdirectories: <indir>/species, <indir>/genes""")
	parser.add_argument('--midas_db', metavar='CHAR', type=str, default=os.environ['MIDAS_DB'] if 'MIDAS_DB' in os.environ else None,
		help="""Path to reference database
By default, the MIDAS_DB environmental variable is used""")
	parser.add_argument('--min_strain_freq', type=float, default=0.10, metavar='FLOAT',
		help="""Minimum relative frequency of strains per species (default=0.10, range=0.0-1.0)
Useful for minimizing spurious, low-frequency predictions""")
	parser.add_argument('--max_similarity', metavar='FLOAT', type=float, default=0.90,
		help="""Maximum gene content between reference strains (default=0.90, range=0.0-1.0)
Used for grouping strains into non-redundant clusters""")
	parser.add_argument('--max_gene_freq', metavar='FLOAT', type=float, default=0.50,
		help="""Maximum frequency of genes across reference strains (default=0.50, range=0.0-1.0)
Used for removing genes present in all reference strains, which are not informative""")
	parser.add_argument('--max_genomes', metavar='INT', type=int, default=300,
		help="""Maximum number of genomes per species (default=300, range=0-inf)
Species with hundreds of genomes currently cause the program to run very slowly""")
	args = vars(parser.parse_args())
	return args

def fetch_gene_sets(midas_db, species_id):
	gene_sets = defaultdict(set)
	inpath = '%s/pan_genomes/%s/gene_info.txt.gz' % (midas_db, species_id)
	for r in csv.DictReader(gzip.open(inpath), delimiter='\t'):
		gene_id = r['centroid_99']
		genome_id = r['genome_id']
		gene_sets[genome_id].add(gene_id)
	return gene_sets

def cluster_strains(midas_db, species_id, cutoff=0.05):
	gene_sets = fetch_gene_sets(midas_db, species_id)
	genome_ids = gene_sets.keys()
	cluster_to_genomes = {}
	genome_to_cluster = {}
	for genome_id in genome_ids:
		cluster_to_genomes[genome_id] = [genome_id]
		genome_to_cluster[genome_id] = genome_id
	for g1, g2 in itertools.combinations(genome_ids, 2):
		d = jaccard(g1, g2, gene_sets)
		if d > cutoff:
			c1, c2 = genome_to_cluster[g1], genome_to_cluster[g2]
			if c1 != c2:
				cluster_to_genomes[c1] += cluster_to_genomes[c2]
				for g in cluster_to_genomes[c2]:
					genome_to_cluster[g] = c1
				del cluster_to_genomes[c2]
	return cluster_to_genomes, genome_to_cluster

def jaccard(g1, g2, gene_sets):
	return 1.0 * len(gene_sets[g1].intersection(gene_sets[g2])) / len(gene_sets[g1].union(gene_sets[g2]))

def read_copy_num(midas_db, species_id, rep_ids):
	copy_num = {}
	inpath = '%s/pan_genomes/%s/gene_info.txt.gz' % (args['midas_db'], species_id)
	reader = csv.DictReader(gzip.open(inpath), delimiter='\t')
	for r in reader:
		genome_id = r['genome_id']
		if genome_id not in rep_ids:
			continue
		if genome_id not in copy_num:
			copy_num[genome_id] = {}
		if r['centroid_99'] not in copy_num[genome_id]:
			copy_num[genome_id][r['centroid_99']] = 0
		copy_num[genome_id][r['centroid_99']] += 1
	return copy_num

def id_genes(copy_num, max_gene_freq):
	gene_prev = {}
	for genome_id in copy_num:
		for gene_id, count in copy_num[genome_id].items():
			if gene_id not in gene_prev:
				gene_prev[gene_id] = 0
			if count > 0:
				gene_prev[gene_id] += 1
	for gene_id, count in gene_prev.items():
		gene_prev[gene_id] = count/len(copy_num)
	gene_ids = sorted([gene_id for gene_id, prev in gene_prev.items() if prev <= max_gene_freq])
	return gene_ids

if __name__ == "__main__":

	args = parse_args()
	
	print "\nchecking arguments"
	if not args['midas_db']:
		error = "\nError: No reference database specified\n"
		error += "Use the flag --db_dir to specify a database,\n"
		error += "Or set the MIDAS_DB environmental variable: export MIDAS_DB=/path/to/midas/db\n"
		sys.exit(error)
	
	print "\nreading genome information"
	genome_info = {}
	reader = csv.DictReader(open('%s/genome_info.txt' % args['midas_db']), delimiter='\t')
	for r in reader:
		genome_info[r['genome_id']] = r

	print "\nreading species information"
	species_info = {}
	reader = csv.DictReader(open('%s/species_info.txt' % args['midas_db']), delimiter='\t')
	for r in reader:
		species_info[r['species_id']] = r

	print "\nreading species abundances"
	species_abundances = {}
	reader = csv.DictReader(open('%s/species/species_profile.txt' % args['indir']), delimiter='\t')
	for r in reader:
		species_abundances[r['species_id']] = r

	species_ids = [_.rstrip() for _ in open('%s/genes/species.txt' % args['indir'])]
	print "\nestimating strain abundances for %s species:" % len(species_ids)
	for species_id in species_ids:
		print "  %s" % species_id

	strain_abundances = []
	for index, species_id in enumerate(species_ids):
		print "\n", "%s. %s" % (index+1, species_id)
		
		if int(species_info[species_id]['count_genomes']) > args['max_genomes']:
			print "skipping species with %s strains" % species_info[species_id]['count_genomes']
			continue
		elif int(species_info[species_id]['count_genomes']) == 1:
			print "only 1 reference strain; nothing to do"
			continue
		
		print "reading read-depth (i.e. coverage) of pangenome genes from MIDAS output"
		coverage = {}
		inpath = '%s/genes/output/%s.genes.gz' % (args['indir'], species_id)
		for r in csv.DictReader(gzip.open(inpath), delimiter='\t'):
			coverage[r['gene_id']] = float(r['copy_number'])
		print "  total genes: %s" % len(coverage)

		print "clustering reference strains at %s%% gene content similarity" % (100*args['max_similarity'])
		cluster_to_genomes, genome_to_cluster = cluster_strains(args['midas_db'], species_id, args['max_similarity'])
		print "  total strains: %s, clustered strains: %s" % (len(genome_to_cluster), len(cluster_to_genomes))

		print "selecting subset of genes found in <=%s%% of reference strains to use for strain estimation" % (100*args['max_gene_freq'])
		copy_num = read_copy_num(args['midas_db'], species_id, cluster_to_genomes)
		gene_ids = id_genes(copy_num, args['max_gene_freq'])
		print "  retained genes: %s" % len(gene_ids)

		print "obtaining initial estimation of strain frequencies"
		genome_ids = sorted(copy_num)
		matrix = np.array([[copy_num[genome_id][gene_id] if gene_id in copy_num[genome_id] else 0 for genome_id in genome_ids] for gene_id in gene_ids])
		vector = np.array([coverage[gene_id] for gene_id in gene_ids])
		coefficients, residual = nnls(A=matrix, b=vector)
		obs_abundances = dict([(g, c/sum(coefficients)) for g, c in zip(genome_ids, coefficients) if c > 0])
		print "  total strains detected: %s" % len(obs_abundances)

		if args['min_strain_freq'] > 0:
			print "re-estimating frequencies after removing strains with <%s%% relative frequency" % (100*args['min_strain_freq'])
			genome_subset = [genome_id for genome_id, obs_abun in obs_abundances.items() if obs_abun >= args['min_strain_freq']]
			matrix = np.array([[copy_num[genome_id][gene_id] if gene_id in copy_num[genome_id] else 0 for genome_id in genome_subset] for gene_id in gene_ids])
			vector = np.array([coverage[gene_id] for gene_id in gene_ids])
			coefficients, residual = nnls(A=matrix, b=vector)
			obs_abundances = dict([(g, c/sum(coefficients)) for g, c in zip(genome_subset, coefficients)])
			print "  total strains detected: %s" % len(obs_abundances)

		tot_abun = float(species_abundances[species_id]['relative_abundance'])
		tot_reads = float(species_abundances[species_id]['count_reads'])
		tot_depth = float(species_abundances[species_id]['coverage'])
		for genome_id, obs_abun in obs_abundances.items():
			row = [species_id, genome_info[genome_id]['genome_name'], genome_id, ','.join(cluster_to_genomes[genome_id]), obs_abun*tot_abun, obs_abun*tot_reads, obs_abun*tot_depth]
			strain_abundances.append(row)

	print "\nwriting results to: %s/strains/strain_freqs.tsv" % args['indir']
	outdir = args['indir']+'/strains'
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	out = open(outdir+'/strain_freqs.tsv', 'w')
	header = ['species_id', 'genome_name', 'genome_id', 'clustered_ids', 'relative_abundance', 'count_reads', 'coverage']
	out.write('\t'.join(header)+'\n')
	for row in strain_abundances:
		out.write('\t'.join([str(_) for _ in row])+'\n')


