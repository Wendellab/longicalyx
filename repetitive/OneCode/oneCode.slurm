#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 168:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 150   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "repeatmasker"   # job name
#SBATCH --mem=300G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-RMredo.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=corrinne@iastate.edu   # email address


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

#mkdir F1.dir

module load repeatmasker/4.0.7-openmpi3-vruhfro
RepeatMasker -pa 150 -s -nolow -lib cottonRB23.04.fa -no_is -frag 500000 -nocut -noisy -dir F1.dir -html -gff F1.fasta

cp F1.dir/F1.fasta.out .

perl build_dictionary.pl --unknown --rm F1.fasta.out > F1.fasta.out.dictionary
perl one_code_to_find_them_all.pl --rm F1.fasta.out --ltr F1.fasta.out.dictionary --unknown 
grep "###" *.copynumber.csv > F1.summary

sed -i 's/F1[.]fasta[.]out_F/\tChr/g' F1.summary
sed -i 's/[.]copynumber[.]csv:\#\#*/\t/g' F1.summary
