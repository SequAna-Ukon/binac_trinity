#!/usr/bin/env nextflow

/* 
The aim of this Nextflow script is 2 fold.
1 - I want to use it to document the full work flow
    of going from multiple pairs of raw seq fastq.gz files
    to the trinity denovo assembly and possibly even metrics of the assmebly.
2 - I want to get an idea of the resource requirments for a given number of reads.

Because of 2 I will work with subsampling to work with 4 or 5 subsets of the
available reads so that we can get an idea of how the processing time scales
with number of reads.

I beleive BL was talking of 20+ sequenced transciptomes.

Much of the code we use here will be taken from the pipeline that I developed with
NY.

I will comment up thoroughly so that this script can be used as a reference.

The original sequencing files that NY had totalled 35G worth of sequencing data.
This is quite a lot of data so we should have plenty to get an idea of run time.

I just counted the reads in 1 of her files and it contains ~45M reads.

So we will work with just two of the pairs of files (so that we show the concatenation)
but not make use of the other files.
*/

// Full path to the directory containing the fastq.gz files
raw_reads_dir = "/home/humebc/projects/20201217_yamada/raw_seq_files"

// This gives us 2 pairs of read files
raw_reads_in_ch = Channel.fromPath("${raw_reads_dir}/*_DC_*_{1,2}.fastq.gz")

// The full path to the SILVA fastq.gz db that has been created by: cat SILVA_138.1_SSUParc_tax_silva_trunc.fasta.gz SILVA_138.1_LSUParc_tax_silva_trunc.fasta.gz > SILVA_138.1_SSUParc.LSUParc.tax_silva_trunc.fasta.gz
silva_fastq_gz = file("/home/humebc/projects/bernard_lepetit/SILVA_db/SILVA_138.1_SSUParc.LSUParc.tax_silva_trunc.fasta.gz")

// The path to the bin directory
remove_unfixable_script = file("/home/humebc/projects/bernard_lepetit/bin/remove_unfixable.py")

// Concatenate the sequencing files together
// N.B. Its important that the files are concatenated together in the same order
// To do this we'll use the find command.
process concatenate_fasta{
    tag "concatenate"
    publishDir "/home/humebc/projects/bernard_lepetit/results/concatenate"
    cpus 1

    input:
    path(fastqs) from raw_reads_in_ch.collect()

    output:
    tuple path("reads.R1.fastq.gz"), path("reads.R2.fastq.gz") into concatenate_out_ch

    script:
    """
    find -L . -type f -name "*_1.fastq.gz" | sort | xargs cat > reads.R1.fastq.gz
    find -L . -type f -name "*_2.fastq.gz" | sort | xargs cat > reads.R2.fastq.gz
    """
}


// subsample the raw reads to a given number of reads
process subsample{
    tag "subsample ${proportion}"
    container 'didillysquat/seqtk_ps:latest'
    publishDir "/home/humebc/projects/bernard_lepetit/results/subsample"
    cpus 1

    input:
    each proportion from Channel.fromList([1000000, 10000000, 50000000])
    tuple path(r1), path(r2) from concatenate_out_ch

    output:
    tuple val(proportion), path("reads.raw.subsample.${proportion}.R1.fastq.gz"), path("reads.raw.subsample.${proportion}.R2.fastq.gz") into subsample_out_ch

    script:
    """
    seqtk sample -s100 $r1 ${proportion} | gzip > reads.raw.subsample.${proportion}.R1.fastq.gz
    seqtk sample -s100 $r2 ${proportion} | gzip > reads.raw.subsample.${proportion}.R2.fastq.gz
    """
}

// QC of the raw reads to remove adapter sequences and remove low quality sequences
process fastp{
    tag "fastp val: ${proportion}"
    container "biocontainers/fastp:v0.20.1_cv1"
    publishDir "/home/humebc/projects/bernard_lepetit/results/fastp"
    cpus 1

    input:
    tuple val(proportion), path(r1), path(r2) from subsample_out_ch

    output:
    file "subsample.${proportion}.fastp.html" into fastp_out_html_ch
    tuple val(proportion), path("reads.subsample.${proportion}.fastp.R1.fastq.gz"), path("reads.subsample.${proportion}.fastp.R2.fastq.gz") into fastp_out_clean_ch

    script:
    """
    fastp -q 20 -i $r1 -I $r2 -o reads.subsample.${proportion}.fastp.R1.fastq.gz -O reads.subsample.${proportion}.fastp.R2.fastq.gz
    mv fastp.html subsample.${proportion}.fastp.html
    """
}

// Now error correction
// Followed by removal of the unfixable_error reads using
// the harvard script that I have modified to work with python3
process rcorrector{
    tag "rcorrector ${proportion}"
    container 'tabotaab/rcorrector:latest'
    publishDir "/home/humebc/projects/bernard_lepetit/results/corrected_reads"
    cpus 28
    
    input:
    tuple val(proportion), path(r1), path(r2) from fastp_out_clean_ch
    path remove_unfixable_script

    output:
    tuple val(proportion), path("reads.subsample.${proportion}.fastp.R1.cor.unfixable_removed.fq.gz"), path("reads.subsample.${proportion}.fastp.R2.cor.unfixable_removed.fq.gz") into rcorrector_out_ch

    script:
    """
    run_rcorrector.pl -1 ${r1} -2 ${r2} -od . -t ${task.cpus}
    python3 $remove_unfixable_script -1 reads.subsample.${proportion}.fastp.R1.cor.fq.gz -2 reads.subsample.${proportion}.fastp.R2.cor.fq.gz -s $proportion
    mv unfixrm_reads.subsample.${proportion}.fastp.R1.cor.fq reads.subsample.${proportion}.fastp.R1.cor.unfixable_removed.fq
    gzip reads.subsample.${proportion}.fastp.R1.cor.unfixable_removed.fq
    mv unfixrm_reads.subsample.${proportion}.fastp.R2.cor.fq reads.subsample.${proportion}.fastp.R2.cor.unfixable_removed.fq
    gzip reads.subsample.${proportion}.fastp.R2.cor.unfixable_removed.fq
    """
}

// Remove rRNA reads
// Discard any read pairs where one or both reads map to the SILVA rRNA database
// Output the metrics of the mapping

// First build the index
process bowtie2_make_silva_index{
    tag "bowtie2_silva_mapping ${proportion}"
    container 'biocontainers/bowtie2:v2.4.1_cv1'
    cpus 28

    input:
    path silva_fastq_gz

    output:
    path("SILVA_138.1_SSUParc.LSUParc.tax_silva_trunc.idx*") into bowtie2_make_silva_index_out_ch

    script:
    """
    bowtie2-build --threads ${task.cpus} --large-index $silva_fastq_gz SILVA_138.1_SSUParc.LSUParc.tax_silva_trunc.idx
    """
}


// Then do the mapping
process bowtie2_silva_mapping{
    tag "bowtie2_silva_mapping ${proportion}"
    container 'biocontainers/bowtie2:v2.4.1_cv1'
    publishDir "/home/humebc/projects/bernard_lepetit/results/silva_mapping_out"
    cpus 28

    input:
    tuple val(proportion), path(r1), path(r2), path(index) from rcorrector_out_ch.combine(bowtie2_make_silva_index_out_ch.map{[it]})
    path silva_fastq_gz

    output:
    tuple val(proportion), path("reads.subsample.${proportion}.fastp.R1.cor.norrna.fq.gz"), path("reads.subsample.${proportion}.fastp.R2.cor.norrna.fq.gz") into bowtie2_silva_mapping_out_ch

    script:
    """
    bowtie2 --nofw --quiet --very-sensitive-local --phred33  -x SILVA_138.1_SSUParc.LSUParc.tax_silva_trunc.idx -1 $r1 -2 $r2 --threads ${task.cpus} \\
    --al-conc-gz reads.subsample.${proportion}.fastp.cor.paired.aligned.fq.gz \\
    --un-conc-gz reads.subsample.${proportion}.fastp.cor.paired.unaligned.fq.gz  --al-gz reads.subsample.${proportion}.fastp.cor.unpaired.aligned.fq.gz \\
    --un-gz reads.subsample.${proportion}.fastp.cor.unpaired.unaligned.fq.gz
    mv reads.subsample.${proportion}.fastp.cor.paired.unaligned.fq.1.gz reads.subsample.${proportion}.fastp.R1.cor.norrna.fq.gz
    mv reads.subsample.${proportion}.fastp.cor.paired.unaligned.fq.2.gz reads.subsample.${proportion}.fastp.R2.cor.norrna.fq.gz
    """
}


// Now we can run an assembly.
process trinity{
    cache 'lenient'
    tag "trinity ${proportion}"
    container 'trinityrnaseq/trinityrnaseq:latest'
    cpus 28
    memory '120 GB'
    publishDir "/home/humebc/projects/bernard_lepetit/results/trinity_assembly"

    input:
    tuple val(proportion), path(r1), path(r2) from bowtie2_silva_mapping_out_ch

    output:
    tuple val(proportion), path("${proportion}.Trinity.fasta") into trinity_assembly_out_ch
    tuple val(proportion), path("${proportion}.gene_trans_map") into trinity_gene_trans_map_out_ch
    path "${proportion}.collectl.dat" into trinity_metrics_out_ch

    script:
    // NB that the output directory for each trinity assembly must have 'trinity' in it.
    """
    Trinity --left $r1 --right $r2 --seqType fq --max_memory 120G --CPU ${task.cpus} \\
    --min_contig_length 250 --output trinity --full_cleanup --SS_lib_type RF --monitoring
    mv trinity.Trinity.fasta ${proportion}.Trinity.fasta
    mv trinity.Trinity.fasta.gene_trans_map ${proportion}.gene_trans_map
    """
}

