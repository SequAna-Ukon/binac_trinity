# down-stream Trinity de-novo assembly on BINAC
## Taxonomy profiling for reads
implement some taxonomic checks to see that we are working with the organisms that we expect An mmseqs taxonomy database has been created as per section "Taxonomy assignment" of:https://mmseqs.com/latest/userguide.pdf
- concatenate reads from Ben's pipeline were used
- concatenate reads were sub-sampled for 10K reads

````bash
seqtk sample -s100 reads.R1.fastq.gz 10000 | gzip > reads_sub10k.R1.fastq.gz
seqtk sample -s100 reads.R2.fastq.gz 10000 | gzip > reads_sub10k.R2.fastq.gz
````
- Taxonomy profiling using  kraken2 against nt database

````bash
#download indexed database from indexs zone
kraken2-build --threads 16 --download-taxonomy --db kraken_nt
kraken2-build --threads 16 --download-library nt --db kraken_nt
kraken2-build --build --threads 16 --db kraken_nt
kraken2 --db kraken_nt --threads 30 --report kraken2_report.tsv --paired reads_sub10k.R1.fastq.gz  reads_sub10k.R2.fastq.gz > /dev/null
sort -k4,4 -k2,2rn kraken2_report.tsv  | grep G | head

#krona

ktUpdateTaxonomy.sh --only-build /share/databases/k2_nt/taxonomy/

ktImportTaxonomy -o krona_report.html -t 5 -m 3 -tax /share/databases/k2_nt/taxonomy/ kraken2_report.tsv 
````


- we follow the downstream analyses in https://github.com/vondrakt/de_novo_RNA-seq_analysis,based on the Trinity assembly done by Ben. 
- i found out that the tuterial read me is defferent that actual pipeline setup, will change it later.
- we will follow the startgy that i develop the command the Bernard run it on BINAC, the tools will installed one by one on mamba enviroment 
- we repaeated the fastp since the individual cleaned reads files are missing from the Ben's nextflow:

````bash
while read i;do fastp -w 26 -q 28 -i raw_reads/raw_seq_data/$i/${i}_1.fq.gz -I raw_reads/raw_seq_data/$i/${i}_2.fq.gz -o trim_reads/${i}_1.trim.fq.gz -O trim_reads/${i}_2.trim.fq.gz && rm trim_reads/fastp.html trim_reads/$i.html;done < samples.txt
````

## Transcriptome Assembly Quality Assessment
-- needs to install star, busco and trinity

- Examine the RNA-Seq read representation of the assembly. Ideally, at least ~80% of your input RNA-Seq reads are represented by your transcriptome assembly. The remaining unassembled reads likely to correspond to lowly expressed transcripts with insufficient coverage to enable assembly, or are low quality or aberrant reads. (will already go through after star alignments)



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
TransDecoder.LongOrfs -t 1.Trinity.fasta
TransDecoder.Predict -t Trinity.fasta
#convert_from_gff_to_gtf (because star and stringtie reuired gtf)
agat_convert_sp_gff2gtf.pl -gff 1.Trinity.fasta.transdecoder.gff3 -o Trinity.fasta.transdecoder_agat.gtf
````

### Transcript mapping
star=2.7.10 was installed

````bash
STAR --runMode genomeGenerate --runThreadN 30 --genomeDir trinity_index --genomeFastaFiles Trinity.fasta --genomeSAindexNbases 10 --sjdbGTFfile Trinity.fasta.transdecoder_agat.gtf

while read i;do STAR --genomeDir Trinity_index --runThreadN 30 --readFilesIn trim_reads/${i}_1.trim.fq.gz trim_reads/${i}_2.trim.fq.gz --readFilesCommand zcat --quantMode GeneCounts --outFileNamePrefix $i --outSAMtype BAM SortedByCoordinate;done < samples.txt 
````

### Transcript Quantification
-- install stringtie2

````bash
while read i;do stringtie $i.sorted.bamAligned.sortedByCoord.out.bam -p 28 -G 1.Trinity.fasta.transdecoder_agat.gtf -e -o $i.gtf -A $i.gene_abundances.tsv;done < samples.txt

wget https://ccb.jhu.edu/software/stringtie/dl/prepDE.py3 

python3 prepDE.py3 -i prepDE_samples.txt
````
- The configuration "prepDE_samples.txt" file must be in format:

sample1[tab]path_to_sample1.gtf

sample2[tab]path_to_sample2.gtf 

sample3[tab]path_to_sample3.gtf

### DE analysis (following Trinotate)
-- install trinity=v2.15.1

- The configuration file "config_DE.txt" for the DE analysis must be in the format:

condition1[tab]sample1

condition1[tab]sample2

condition1[tab]sample3

condition2[tab]sample4

condition2[tab]sample5

condition2[tab]sample6

- convert the matrix from csv to matrix
````bash
sed 's/,/\t/g' gene_count_matrix.csv > gene_count.matrix
OR
sed 's/,/\t/g' transcript_count_matrix.csv > transcript_count.matrix
````
#### Tips before run the analysis:
  
-	Create and activate new environment
  
-	Install latest R:

mamba install r-base=4.2.2

mamba install -c conda-forge r-curl

- Start R:
R

-	Install BiocManager:
  
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
 	
install.packages("survival")

BiocManager::install("DESeq2")

BiocManager::install("edgeR")

BiocManager::install("cluster")

BiocManager::install("qvalue")

BiocManager::install("fastcluster")

BiocManager::install("Biobase")

-	Validate that DESeq2 and edgeR are working:

library(DESeq2)

library(edgeR)

- Exit R:
Ctrl+D

- download the latest Trinity from GitHub  

git clone https://github.com/trinityrnaseq/trinityrnaseq.git

- export the path for the Trinity scripts

export PATH=$PATH:trinityrnaseq/Analysis/DifferentialExpression/

#### edgeR


````bash
run_DE_analysis.pl --matrix transcript_count.matrix --samples_file config_DE.txt --reference_sample T0 --method edgeR --output edgeR_transcript
````

#### DESeq2

````bash
run_DE_analysis.pl --matrix transcript_count.matrix  --samples_file config_DE.txt --reference_sample T0 --method DESeq2 --output DESeq2_transcript
````
#### subsetting DE

- for edgeR and DESeq2 and with 4 fold of change and Pvalue 0.05 

````bash
cd DESeq2_transcript
analyze_diff_expr.pl --matrix ../transcript_count.matrix --samples ../config_DE.txt -P 0.05 -C 2
define_clusters_by_cutting_tree.pl -R  diffExpr.P0.05_C2.matrix.RData --Ptree 60
````


### Trinotate Functional Annotation

- notes from Ben about singularity in BINAC
To run singularity on BINAC, you'll need to install it in your Home File System. Then you'll need to set the SINGULARITY_CACHE_DIR environmental variable, for example: export SINGULARITY_CACHEDIR=/beegfs/work/kn_pop528802/.singularity It is a good idea to put this line in your .bashrc file so that you don't have to remember to set it each time. N.B. your SINGULARITY_CACHE_DIR cannot be in your Home File Sytem but must be in your Work File System.


#will be using Trinotate through docker

````bash
docker pull trinityrnaseq/trinotate

docker run --rm -it -v `pwd`:/data -v /tmp:/tmp -e TRINOTATE_HOME=/usr/local/src/Trinotate trinityrnaseq/trinotate bash

export PATH=$PATH:/usr/local/src/Trinotate/
````
#TMHMM and  signalp installed sepra.

- download TMHMM

https://services.healthtech.dtu.dk/services/TMHMM-2.0/
````bash

tar -xzf tmhmm-2.0c.Linux.tar.gz

which perl

add perl path in tmhmm-2.0c/bin/tmhmmformat.pl and tmhmm-2.0c/bin/tmhmm

export PATH=$PATH:tmhmm-2.0c/bin/

tmhmm --short 1.Trinity.fasta.transdecoder.pep > tmhmm.v2.out

````
- download signalp-6
  
https://services.healthtech.dtu.dk/service.php?SignalP-6.0.

````bash

tar -xzf signalp-6.0h.fast.tar.gz

pip install signalp-6-package/

SIGNALP_DIR=$(python3 -c "import signalp; import os; print(os.path.dirname(signalp.__file__))" )

cp -r signalp-6-package/models/* $SIGNALP_DIR/model_weights/

signalp6 --fastafile 1.Trinity.fasta.transdecoder.pep --output_dir sigP6outdir --format none --organism euk --mode fast

````
- Trinotate commands
````bash

export PATH=/data/eggnog-mapper:/data/eggnog-mapper/eggnogmapper/bin:"$PATH"
export PATH=$PATH:/data/eggnog-mapper

Trinotate --create --db cap_Trinotate.sqlite --trinotate_data_dir cap_Trinotate --use_diamond

Trinotate --db cap_Trinotate.sqlite --init --gene_trans_map 1.gene_trans_map --transcript_fasta 1.Trinity.fasta --transdecoder_pep 1.Trinity.fasta.transdecoder.pep 

Trinotate --db cap_Trinotate.sqlite --LOAD_tmhmmv2 tmhmm.v2.out

Trinotate --db cap_Trinotate.sqlite --LOAD_signalp sigP6outdir/output.gff3


Trinotate --db cap_Trinotate.sqlite --CPU 50 --transcript_fasta 1.Trinity.fasta --transdecoder_pep 1.Trinity.fasta.transdecoder.pep --trinotate_data_dir cap_Trinotate --run "swissprot_blastp swissprot_blastx pfam infernal EggnogMapper" --use_diamond

Trinotate --db cap_Trinotate.sqlite --report > cap_Trinotate.xls

````


- covert annotation file to html

````bash
install.packages('remotes') 
remotes::install_github('rstudio/DT')
library(DT)
data <- as.matrix(read.table("cap_Trinotate.xls", sep = '\t', header = TRUE, fill = TRUE))
Dqua_DT <- datatable(data)
saveWidget(Dqua_DT, "cap_Trinotate.html", selfcontained = TRUE, libdir = NULL, title = "Captiva.annotations", knitrOptions = list())
````
- GSEA using goatools
https://github.com/tanghaibao/goatools/tree/main

wget http://current.geneontology.org/ontology/go-basic.obo
wget http://current.geneontology.org/ontology/subsets/goslim_generic.obo

pip install goatools

- preparing of assioation/gene2go file "trans_GO_egg.txt"
  cut -f1,10 cap_Trinotate.xls| sed 's/,/;/g' > trans_GO_egg.txt
- preparing of study file 
- just list of the target DE subsets
- preparing of population file "trans_all.ids"
  grep '>' 1.Trinity.fasta.transdecoder.pep |sed 's/>//g' > trans_all.ids
  
- run find_enrichment.py script on each subset
find_enrichment.py T0_vs_T180_B_Down.ls trans_all.ids  ../trans_GO_egg.txt  --outfile T0_vs_T180_B.depl.xls

- ".depl" for depleted and ".enr" for enriched.
- you can filter on the significance of (e)nrichment or (p)urification. it can report various multiple testing corrected p-values as well as the false discovery rate. The e in the "Enrichment" column means "enriched" - the concentration of GO term in the study group is significantly higher than those in the population. The "p" stands for "purified" - significantly lower concentration of the GO term in the study group than in the population.
