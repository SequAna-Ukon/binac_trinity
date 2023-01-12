#PBS -l nodes=1:ppn=5
#PBS -l walltime=00:05:00
#PBS -l mem=5gb
#PBS -S /bin/bash
#PBS -N fastp_with_conda_test
#PBS -j oe
#PBS -o fastp_with_conda.log
cd /beegfs/work/kn_pop528802/testing/fastp_with_conda 
source /home/kn/kn_kn/kn_pop528802/local_software/miniconda3/etc/profile.d/conda.sh
conda activate fastp_env
fastp -q 20 -i reads.raw.subsample.10000.R1.fastq.gz -I reads.raw.subsample.10000.R2.fastq.gz -o reads.raw.subsample.10000.R1.fastp.fastq.gz -O reads.raw.subsample.10000.R2.fastp.fastq.gz 
