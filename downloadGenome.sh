#!/bin/bash 
#SBATCH -p defq
#SBATCH --mem-per-cpu=4G
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --output OUT/%J
#SBATCH --error ERR/%J 
#SBATCH --job-name=genome
#SBATCH --account=sbi9srj

mkdir -p resources/humanGRCh38
cd resources/humanGRCh38

#Download FASTA
curl -sLO https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz

#Download GFF
curl -sLO https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz

gunzip *
