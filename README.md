# Snakemake workflow: 2023_ecoli_protection

## Quick start

Create a symbolic links to a copy of the Uniref50 protein fasta file, changing the following command
to the actual PATH to that file on your system:

    ln -s /storage/datasets/uniprot/uniref50.fasta data/uniref50.fasta

Tip: the file [can be downloaded from this page](https://www.uniprot.org/help/downloads) ([direct link to the file](https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz))

Also create a symbolic link to the directory in which the databases for eggnog-mapper have been placed:

    ln -s /fast-storage/miniconda3/envs/eggnog-mapper/lib/python3.9/site-packages/data/ data/eggnog-mapper

Note: the above path would likely be different in your system. The best way to get these files is to install
`eggnog-mapper` using a `conda` environment and then use the `download_eggnog_data.py` command.

Place your genome assemblies in the `data/gffs` and `data/fastas` directories. The files should be named
`SAMPLE.gff` and `SAMPLE.fasta`, respectively, and the sample names should match those in the phenotype file
(i.e. `data/data.tsv`).

Create and activate a `conda` environment to run the bootstrapping script and the pipeline (named `microGWAS`, can be skipped if it's already present):

    conda env create -f environment.yml
    conda activate microGWAS

Then run the bootstrapping script to populate the input files for the pipeline and download the reference genomes
used for annotation of hits and the rare variants analyses.

    bash bootstrap.sh Escherichia coli IAI39 GCF_000013305.1,GCF_000007445.1,GCF_000026305.1,GCF_000026265.1,GCF_000026345.1,GCF_000005845.2,GCF_000026325.1,GCF_000013265.1 

You are now ready to run the full pipeline! The following example runs all the analyses using 24 cores and `mamba` as the conda backend
to install each environment:

    snakemake -p annotate_summary panfeed map_back heritability enrichment qq_plots find_amr_vag tree --cores 24 --use-conda --conda-frontend mamba
    
## Author

Marco Galardini (marco.galardini@twincore.de)
