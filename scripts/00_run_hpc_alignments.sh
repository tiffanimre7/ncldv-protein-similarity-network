#!/bin/bash
#SBATCH --job-name=ncldv_psn
#SBATCH --cpus-per-task=8
#SBATCH --mem=30000M
#SBATCH --time=06:00:00
#SBATCH --output=ncldv_psn_%j.out
#SBATCH --error=ncldv_psn_%j.err

set -e

echo "Job started: $(date)"
echo "Working directory before temp setup: $(pwd)"
echo "Running on host: $(hostname)"
echo "SLURM job ID: $SLURM_JOB_ID"
echo "Using 8 threads, matching the SLURM CPU request."

echo "Creating temporary working directory exactly as in the practical..."
mkdir -p /tmp/$USER
cd /tmp/$USER
echo "Now working in temporary directory: $(pwd)"

echo "Copying the FULL Nucleocytoviricota FASTA file, not the partial demo file..."
cp /lisc/data/scratch/course/2026s053521/01psn/Nucleocytoviricota.faa ./ncldv.faa

echo "Checking number of protein sequences in the copied FASTA..."
grep -c "^>" ncldv.faa

echo "Creating BLAST protein database..."
module load BLAST+
module list BLAST+
makeblastdb -dbtype prot -in ncldv.faa -parse_seqids -blastdb_version 4
module unload BLAST+

echo "Running SWIPE all-vs-all local alignment..."
module load SWIPE
module list SWIPE
/usr/bin/time -v -o ~/assignment_sequence_alignments/ncldv.swipe.time.txt \
    swipe \
    --query ncldv.faa \
    --db ncldv.faa \
    --out ncldv.swipe.m8 \
    --outfmt 8 \
    --matrix BLOSUM45 \
    --gapopen 11 \
    --gapextend 3 \
    --evalue 1
module unload SWIPE

echo "Running BLASTP all-vs-all local alignment..."
module load BLAST+
module list BLAST+
/usr/bin/time -v -o ~/assignment_sequence_alignments/ncldv.blastp.time.txt \
    blastp \
    -query ncldv.faa \
    -db ncldv.faa \
    -out ncldv.blastp.m8 \
    -outfmt 6 \
    -matrix BLOSUM45 \
    -gapopen 11 \
    -gapextend 3 \
    -evalue 1 \
    -num_threads 8
module unload BLAST+

echo "Removing self hits from SWIPE and BLASTP outputs..."
cat ncldv.swipe.m8 | sed 's/lcl|//' | awk '$1!=$2' > ncldv.swipe.wo_self.m8
cat ncldv.blastp.m8 | sed 's/lcl|//' | awk '$1!=$2' > ncldv.blastp.wo_self.m8

echo "Running MMseqs2 all-vs-all local alignment..."
module load MMseqs2
module list MMseqs2
mmseqs createdb ncldv.faa ncldv
/usr/bin/time -v -o ~/assignment_sequence_alignments/ncldv.mmseqs.time.txt \
    mmseqs search ncldv ncldv result /tmp/$USER --threads 8
mmseqs convertalis ncldv ncldv result ncldv.mmseqs.m8
cat ncldv.mmseqs.m8 | sed 's/lcl|//' | awk '$1!=$2' > ncldv.mmseqs.wo_self.m8
module unload MMseqs2

echo "Moving final self-hit-filtered m8 files to ~/assignment_sequence_alignments..."
mv *.wo_self.m8 ~/assignment_sequence_alignments/

echo "Counting final output lines..."
wc -l ~/assignment_sequence_alignments/*.wo_self.m8 > ~/assignment_sequence_alignments/ncldv_m8_line_counts.txt
cat ~/assignment_sequence_alignments/ncldv_m8_line_counts.txt

echo "Runtime files saved:"
ls -lh ~/assignment_sequence_alignments/*.time.txt

echo "Cleaning up temporary directory exactly as in the practical..."
cd
rm -rf /tmp/$USER

echo "Job finished successfully: $(date)"
