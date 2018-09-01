# MIDAS-strains
Estimate strains from reads mapped to pan-genomes from the MIDAS database


## Usage

Download the new package using github; one new python library `scipy` is required if not already installed:
`git clone https://github.com/snayfach/MIDAS-strains`

Note that the files `reads_1.fa` and and `reads_2.fa` are simulated FASTA files containing 20 strains from 5 different human gut species at equal relative abundances (5%)

MIDAS should already be installed; run MIDAS as usual for the species module to estimate species abudnances:
```
export PYTHONPATH=$PYTHONPATH:/path/to/MIDAS
export PATH=$PATH:/path/to/MIDAS/scripts
export MIDAS_DB=/path/to/midas_db_v1.2

run_midas.py species midas_out -1 reads_1.fa -2 reads_2.fa -t 10
```

Now we use MIDAS to map reads against the pangenomes. Here we can specify the 5 species that we know are in the community to force MIDAS to map reads against these pangenomes regardless of their coverage in the test dataset. Also, I found slighly better performance in the test dataset when increasing the mapping % identity to 99%:
`run_midas.py genes midas_out -1 MIDAS-strains/reads_1.fa -2MIDAS-strains/reads_2.fa -t 10 --mapid 99 --species_id Bacteroides_fragilis_54507,Bacteroides_vulgatus_57955,Bifidobacterium_longum_57796,Clostridium_perfringens_56840,Parabacteroides_distasonis_56985
`

Now, use the script `MIDAS-strains/midas_strains.py` to estimate strain frequencies from the MIDAS output:
`MIDAS-strains/midas_strains.py --indir midas_out`

You should see some output that looks something like this:
```
estimating strain abundances for 5 species:
  Clostridium_perfringens_56840
  Bacteroides_vulgatus_57955
  Bifidobacterium_longum_57796
  Bacteroides_fragilis_54507
  Parabacteroides_distasonis_56985

1. Clostridium_perfringens_56840
reading read-depth (i.e. coverage) of pangenome genes from MIDAS output
  total genes: 17315
clustering reference strains at 90.0% gene content similarity
  total strains: 12, clustered strains: 12
selecting subset of genes found in <=50.0% of reference strains to use for strain estimation
  retained genes: 17048
obtaining initial estimation of strain frequencies
  total strains detected: 9
re-estimating frequencies after removing strains with <10.0% relative frequency
  total strains detected: 4
...
```

With an output file like this:
```
species_id      genome_name     genome_id       clustered_ids   relative_abundance      count_reads     coverage
Clostridium_perfringens_56840   Clostridium perfringens D str. JGS1721  488537.5        488537.5        0.0609780725004 22.5772616892   0.230069514784
Clostridium_perfringens_56840   Clostridium perfringens E str. JGS1987  451755.5        451755.5        0.0638524275139 23.6414977772   0.240914420766
Clostridium_perfringens_56840   Clostridium perfringens ATCC 13124      195103.10       195103.10       0.0612481361613 22.677253337    0.231088461643
Clostridium_perfringens_56840   Clostridium perfringens NCTC 8239       451757.5        451757.5        0.0569989611244 21.1039871966   0.215056376683
```

Finally, you can use the -h flag to see a few options, but I've already set good defualt values, so this can be mostly ignored:
```
./midas_strains.py -h
description: use reads mapped to pangenomes to estimate strain frequencies

usage: midas_strains.py --indir <indir> [options]

optional arguments:
  -h, --help            show this help message and exit
  --indir CHAR          Path to input directory with MIDAS output for one sample
                        Should contain subdirectories: <indir>/species, <indir>/genes
  --midas_db CHAR       Path to reference database
                        By default, the MIDAS_DB environmental variable is used
  --min_strain_freq FLOAT
                        Minimum relative frequency of strains per species (default=0.10, range=0.0-1.0)
                        Useful for minimizing spurious, low-frequency predictions
  --max_similarity FLOAT
                        Maximum gene content between reference strains (default=0.90, range=0.0-1.0)
                        Used for grouping strains into non-redundant clusters
  --max_gene_freq FLOAT
                        Maximum frequency of genes across reference strains (default=0.50, range=0.0-1.0)
                        Used for removing genes present in all reference strains, which are not informative
  --max_genomes INT     Maximum number of genomes per species (default=300, range=0-inf)
                        Species with hundreds of genomes currently cause the program to run very slowly
```

