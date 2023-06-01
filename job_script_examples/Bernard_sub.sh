
#PBS -l nodes=1:ppn=14
#PBS -l walltime=4:00:00
#PBS -l mem=100gb
#PBS -S /bin/bash
#PBS -N Fastp_individuals
#PBS -j oe
#PBS -o Fastp_individuals.log
cd /beegfs/work/kn_pop243393/raw_reads/raw_seq_data
source /home/kn/kn_kn/kn_pop243393/miniconda3/etc/profile.d/conda.sh
conda activate DE
