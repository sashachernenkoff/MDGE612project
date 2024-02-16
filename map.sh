#!/bin/bash
# Slurm Script Input Variables
#SBATCH --job-name=mdge612
#SBATCH --chdir=/work/long_lab/sasha/mdge612
#SBATCH --error=mdge612_%A-%a.error
#SBATCH --output=mdge612_%A-%a.out
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --partition=cpu2022,theia,cpu2021,cpu2019

# Print
echo "=============="
echo "Running Jawamix5"
echo "=============="
echo ""

# Command
start=$SECONDS

# Convert phenotype file to tab delimited format
awk -F, 'BEGIN {OFS="\t"} {print $1, $2}' filtered.FT10.txt > filtered.FT10.tsv

# Convert genotype file to hdf5
java -Xmx2g -jar Jawamix5/jawamix5.jar import -ig filtered.coded_call_method_54.tair9.FT10.csv -o geno.hdf5

# Create kinship matrix
java -Xmx2g -jar Jawamix5/jawamix5.jar kinship -ig geno.hdf5 -o kinship -scale 1 -m RRM

# Run linear regression for phenotype
java -Xmx2g -jar Jawamix5/jawamix5.jar lm -ig geno.hdf5 -ip filtered.FT10.tsv -o lm_out
java -Xmx2g -jar Jawamix5/jawamix5.jar emmax_stepwise -ig geno.hdf5 -ip filtered.FT10.tsv -o emmax_out -ik kinship.RRM

end=$SECONDS
echo "duration: $((end-start)) seconds."