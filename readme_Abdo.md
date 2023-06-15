# down-stream Trinity de-novo assembly on BINAC

- we follow the downstream analyses in https://github.com/vondrakt/de_novo_RNA-seq_analysis,based on the Trinity assembly done by Ben. 
- i found out that the tuterial read me is defferent that actual pipeline setup, will change it later.
- we will follow the startgy that i develop the command the Bernard run it on BINAC, the tools will installed one by one on mamba enviroment 
- we repaeated the fastp since the individual cleaned reads files are missing from the Ben's nextflow:

````bash
while read i;do fastp -w 26 -q 28 -i raw_reads/raw_seq_data/$i/${i}_1.fq.gz -I raw_reads/raw_seq_data/$i/${i}_2.fq.gz -o trim_reads/${i}_1.trim.fq.gz -O trim_reads/${i}_2.trim.fq.gz && rm trim_reads/fastp.html trim_reads/$i.html;done < samples.txt
````

## Transcriptome Assembly Quality Assessment
-- needs to install star, busco and trinity

- Examine the RNA-Seq read representation of the assembly. Ideally, at least ~80% of your input RNA-Seq reads are represented by your transcriptome assembly. The remaining unassembled reads likely corresponds to lowly expressed transcripts with insufficient coverage to enable assembly, or are low quality or aberrant reads. (will already go through after star alignments)



- Compute the E90N50 transcript contig length - the contig N50 value based on the set of transcripts representing 90% of the expression data. for more details (https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats)

````bash
TrinityStats.pl  1.Trinity.fasta
or 
#after the quantification
contig_ExN50_statistic.pl transcripts.TMM.EXPR.matrix 1.Trinity.fasta transcript | tee ExN50.transcript.stats
````

- Use BUSCO to explore completeness according to conserved ortholog content. 

#assuming we working with diatoms

````bash
busco  -i 1.Trinity.fasta -m transcriptome  --cpu 20 -l stramenopiles_odb10 -o busco_stramenopiles
and 
busco  -i 1.Trinity.fasta -m transcriptome --cpu 20 -l eukaryota_odb10 -o busco_eukaryota
````

#we can run it on the Predict coding sequences as well (better)

````bash
busco  -i 1.Trinity.fasta.transdecoder.pep -m protein  --cpu 20 -l stramenopiles_odb10 -o busco_pp_stramenopiles
and 
busco  -i 1.Trinity.fasta.transdecoder.pep -m protein --cpu 20 -l eukaryota_odb10 -o busco_pp_eukaryota
````
#plotting busco report

mkdir busco_plot
cp busco_pp_eukaryota/short*.txt busco_eukaryota/short*.txt busco_stramenopiles/short*.txt busco_pp_stramenopiles/short*.txt
generate_plot.py -wd busco_plot/

## Downstream Analyses

### Predict coding sequences (using TransDecoder)
-- install TransDecoder and agat

````bash
TransDecoder.LongOrfs -t Trinity.fasta
TransDecoder.Predict -t Trinity.fasta
#convert_from_gff_to_gtf (because star and stringtie reuired gtf)
agat_convert_sp_gff2gtf.pl -gff 1.Trinity.fasta.transdecoder.gff3 -o Trinity.fasta.transdecoder_agat.gtf
````

### Transcript mapping
star=2.7.10 was installed

````bash
STAR --runMode genomeGenerate --runThreadN 30 --genomeDir trinity_index --genomeFastaFiles Trinity.fasta --genomeSAindexNbases 10 --sjdbGTFfile Trinity.fasta.transdecoder_agat.gtf

while read i;do STAR --genomeDir Trinity_index --runThreadN 30 --readFilesIn trim_reads/${i}_1.trim.fq.gz trim_reads/${i}_2.trim.fq.gz --readFilesCommand zcat --quantMode GeneCounts --outFileNamePrefix $i.sorted.bam --outSAMtype BAM SortedByCoordinate;done < samples.txt 
````

### Transcript Quantification
-- install stringtie2

````bash
while read i;do stringtie $i.sorted.bam -p 30 -G Trinity.fasta.transdecoder_agat.gtf -e -o $i.gtf -A $i.gene_abundances.tsv;done < samples.txt

wget https://ccb.jhu.edu/software/stringtie/dl/prepDE.py3 

python3 prepDE.py3 -i prepDE_samples.txt
````
- The configuration "prepDE_samples.txt" file must be in format:

sample1[tab]path_to_sample1.gtf

sample2[tab]path_to_sample2.gtf 

sample3[tab]path_to_sample3.gtf

### DE analysis (following Trinotate)
-- install Trinotate

- The configuration file "config_DE.txt" for the DE analysis must be in format

condition1[tab]sample1

condition1[tab]sample2

condition1[tab]sample3

condition2[tab]sample4

condition2[tab]sample5

condition2[tab]sample6


#### edgeR

````bash
run_DE_analysis.pl --matrix gene_count_matrix.tsv or transcript_count_matrix.tsv --samples_file config_DE.txt --reference_sample condition? --method edgeR --output edgeR_genes/trans
````

#### DESeq2

````bash
run_DE_analysis.pl --matrix gene_count_matrix.tsv or transcript_count_matrix.tsv --samples_file config_DE.txt --reference_sample condition? --method DESeq2 --output DESeq2_genes/trans
````
#### subsetting DE

- for edgeR and DESeq2 and with 4 fold of chaneg and Pvalue 0.05 

````bash
analyze_diff_expr.pl --matrix ../gene_count_matrix.tsv --samples config_DE.txt -P 0.05 -C 2
````


### Trinotate Functional Annotation

#will using Trinotate through singularity

````bash
wget https://data.broadinstitute.org/Trinity/TRINOTATE_SINGULARITY/trinotate.v4.0.1.simg
singularity shell -e trinotate.v4.0.1.simg  
export TRINOTATE_HOME=/usr/local/src/Trinotate
$TRINOTATE_HOME/Trinotate
````
#TMHMM and  signalp installed sepra.
pip3 install pybiolib
biolib run --local DTU/DeepTMHMM --fasta input.fasta

https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md

