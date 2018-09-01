# MIDAS-strains
Estimate strains from reads mapped to pan-genomes from the MIDAS database


## Usage

First run MIDAS as normal for the species and genes module:
```
export PYTHONPATH=$PYTHONPATH:/path/to/MIDAS
export PATH=$PATH:/path/to/MIDAS/scripts
export MIDAS_DB=/path/to/midas_db_v1.2

run_midas.py species midas_out -1 reads_1.fa -2 reads_2.fa -t 10
```

I found slighly better performance when increasing the mapping % identity to 99%:
`run_midas.py genes midas_out -1 reads_1.fa -2 reads_2.fa -t 10 --mapid 99`

Download the new package using github:
`git clone https://github.com/snayfach/MIDAS-strains`

One new python library `scipy` is required if not already installed.

Now, use the script `midas_strains.py` to estimate strain frequencies from the MIDAS output:
`/path/to/midas_strains.py --indir midas_out`

You can use the -h flag to see a few options, but I've already set good defualt values, so this can be mostly ignored:
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

