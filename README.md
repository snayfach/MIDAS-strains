# MIDAS-strains
Estimate strains from reads mapped to pan-genomes from the MIDAS database


## Usage

Download the new package using github; one new python library `scipy` is required if not already installed:  
`git clone https://github.com/snayfach/MIDAS-strains`

Note that the files `reads_1.fa` and and `reads_2.fa` are simulated FASTA files containing 20 strains from 5 different human gut species at equal relative abundances (5%). The `abundances.tsv` file lists the exact strains.

MIDAS should already be installed; run MIDAS as usual for the species module to estimate species abudnances:  
```
export PYTHONPATH=$PYTHONPATH:/path/to/MIDAS
export PATH=$PATH:/path/to/MIDAS/scripts
export MIDAS_DB=/path/to/midas_db_v1.2

run_midas.py species midas_out -1 reads_1.fa -2 reads_2.fa -t 10
```

Now we use MIDAS to map reads against the pangenomes. 
Here we can specify the 5 species that we know are in the community to force MIDAS to map reads against these pangenomes regardless of their coverage in the test dataset. 
Also, I found slighly better performance in the test dataset when increasing the mapping % identity to 99%:  
```
run_midas.py genes midas_out -1 MIDAS-strains/reads_1.fa -2MIDAS-strains/reads_2.fa -t 10 --mapid 99 --species_id Bacteroides_fragilis_54507,Bacteroides_vulgatus_57955,Bifidobacterium_longum_57796,Clostridium_perfringens_56840,Parabacteroides_distasonis_56985
```

Now, use the new script to estimate strain frequencies from the MIDAS output:  
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
Bacteroides_fragilis_54507      Bacteroides fragilis str. S6R6  1339321.3       1339321.3,1339320.3,1339322.3,1339319.3,1339317.3       0.0514508615857 19.3461962084
   0.194123465614
Bacteroides_fragilis_54507      Bacteroides fragilis str. S24L26        1339324.3       1339324.3,1339325.3,1339323.3   0.0480365693318 18.0623777101   0.181241383088
Bacteroides_fragilis_54507      Bacteroides sp. 2_1_56FAA       665938.3        665938.3        0.0521754004368 19.6186322832   0.196857141754
Bacteroides_fragilis_54507      Bacteroides fragilis str. 1007-1-F #7   1339292.3       1339292.3,1339293.3,1339340.3,1339295.3,1339294.3,1339339.3,1339338.3,1339337.3 0.0504578045777 18.9727937982   0.190376673782
```

<b>Please notice the column `clustered_ids`. This is important!!</b> This represents redundant strains that are very, very similar. For sample, the B. fragilis strain listed above (genome_name = Bacteroides fragilis str. S6R6 and genome_id = 1339321.3) is estimated to be present at 0.051 (5.1%) relative abundance in the test dataset. But this exact strain is not present in the test dataset (see abundances.tsv). However the strain with genome_id = 1339317.3 is listed in abundances.tsv and is clustered together with 1339321.3. So the output we get is consistent with what we expect.  

For a full mapping of genome ids to genome names, see the MIDAS file: `/path/to/midas_db_v1.2/genome_info.txt`


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

