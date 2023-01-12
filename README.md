# Running a Trinity assembly on BINAC

I have written this guide as a supplement to the guides that are already available for BINAC and to add specific details with regards to running Trinity.

In particular, please see this [quick start guide](https://wiki.bwhpc.de/e/BinAC/Quickstart_Guide) for BINAC. I recommend reading through it.

Please also see the [Trinity documentation](https://github.com/trinityrnaseq/trinityrnaseq/wiki). There are also several tutorials on the Trinity Wiki.

Finally, there is [this](https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html) excellent guide on denovo assemblies from Harvard. The pipelines I have developed to do denovo assemblies, follow the steps described in the guide.

## Connecting to BINAC
You will need to ssh into BINAC with your username and its IP.

For me, this is

`ssh kn_pop528802@login01.binac.uni-tuebingen.de`

It will ask you for your one time password. For me this is a second factor authentication code that I get from the Authenticator app on my phone. Finally it asks for your password.

Login nodes are [here](https://wiki.bwhpc.de/e/BinAC/Login)

## Home File System
Once you are connected you will be in your home directory of the 'Home File System':

```
(base) [kn_pop528802@login01 ~]$ pwd
/home/kn/kn_kn/kn_pop528802
```

From the quickstart guide:
> Home directories are meant for permanent file storage of files that are keep being used like source codes, configuration files, executable programs, conda environments, etc. It is backuped daily and has a quota. If that quota is reached, you will usually experience problems when working on the cluster.


In other words, this is NOT where you store data. This is where you install programs that you want to run.

## Work File System

For me, my work file system directory is:
```
(base) [kn_pop528802@login01 kn_pop528802]$ pwd
/beegfs/work/kn_pop528802
```

From the quickstart guide:
> Use the work file system and not your home directory for your calculations and data.

In other words, this is where you put your data.

## .bashrc
.bashrc is read and sourced every time you log in.

If it doesn't already exist in your home directory, you can create a `.bashrc` file. This behaves the same way as a standard linux .bashrc file. For example you can use this to set your PATH variable. One thing you may want to add to this file is the line:

`module load system/moab`

This loads the system/moab module which is required for some utility programs like checkjob, showq, and showstart (see [here](https://wiki.bwhpc.de/e/BinAC/Specific_Batch_Features)).

To find out more about the .bashrc file, google 'linux .bashrc'.

## Submitting jobs

Jobs are submitted using the qsub command. See the quick start user guide for more.

## Creating job scripts
To run your jobs you will create .sh job scripts that you will then submit using qsub.

You will need to create the job scripts using a text editor. One popular one that should already be availale on BINAC is [nano](https://linuxize.com/post/how-to-use-nano-text-editor/).

See the quick start guide for the basics of what a job script looks like.

I have also created you two example job scripts for running fastp on a set of sequncing files.
One uses [conda](job_script_examples/fastp_with_conda/fastp_with_conda.sh), the other uses [singularity](job_script_examples/fastp_with_singularity/fastp_with_singularity.sh).

## Conda

One convenient way to run programs on BINAC is to install them in a [conda environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands) and then load that environment in your submitted script. E.g. Trinity can be installed with: `conda install -c bioconda trinity` once [conda](https://docs.conda.io/en/latest/miniconda.html) has been installed.

Install miniconda not anaconda.

use 'wget' to download files i.e. 'wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh'

This is probably the easiest way for you to run Trimmomatic and any other programs you want to run, but let's discuss.

Please see the specific instructions for working with conda in your run scripts in the [BINAC documentation](https://wiki.bwhpc.de/e/BinAC/Quickstart_Guide#Queue.2FJob_Basics) (section 3.2) for further info.

[conda environment documentation](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)

## Singularity
Another great option for running programs in BINAC is [Singualrity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html).

Warning: Singularity is a large program with many use options and great capabilities.

Singularity shares many similairties with [Docker](https://www.docker.com/).

Both are means of generating and running [containers](https://www.docker.com/resources/what-container/).

Docker cannot be run on BINAC due to security reasons (it grants root access).

Singularity can be run. Singularity can run Docker '[images](https://jfrog.com/knowledge-base/a-beginners-guide-to-understanding-and-building-docker-images/)'.

A common way to run Singularity is to use it to create a container from a Docker image. There are many Docker images available for many common bioinformatic programs.

You can find them on [DockerHub](https://hub.docker.com/). They are free and you don't need an account to search and download them.

To run singularity on BINAC, you'll need to install it in your Home File System. Then you'll need to set the SINGULARITY_CACHE_DIR environmental variable, for example:
`export SINGULARITY_CACHEDIR=/beegfs/work/kn_pop528802/.singularity`
It is a good idea to put this line in your .bashrc file so that you don't have to remember to set it each time. *N.B.* your SINGULARITY_CACHE_DIR cannot be in your Home File Sytem
but must be in your Work File System.

There are two common ways to run Singularity.

The first is to create an interactive container. You can think about this like booting up a new computer with the programs you want already installed.

E.g. to start a container containing fastp, go to DockerHub and find the desired container then: `singularity shell docker://biocontainers/fastp:v0.20.1_cv1`

You can then work in this terminal using bash as normal.

The second way to run Singularity is to supply a command that you want to run. E.g. for fastp on two files:

`singularity exec docker://biocontainers/fastp:v0.20.1_cv1 fastp -q 20 -i reads.raw.subsample.10000.R1.fastq.gz -I reads.raw.subsample.10000.R2.fastq.gz -o reads.raw.subsample.10000.R1.fastp.fastq.gz -O reads.raw.subsample.10000.R1..fastp.fastq.gz`

This second way is what you will most commonly use on BIANC. You can look at the [example job script](job_script_examples/fastp_with_singularity/fastp_with_singularity.sh) I made for an example.

## Nextflow

[Nextflow](https://www.nextflow.io/) is a program for writing pipelines. I have written a pipeline for doing a denovo assmebly reference [binac_trinity.nf](binac_trinity.nf).

The purpose of this pipeline was 2 fold.

1. I wanted to run this pipeline with varying degrees of subsampling to get an idea of how many resources each of the steps of the pipeline uses.

2. I wanted you to have it as a resource to show the processes of creating a denovo assemply, which containers I use, and how the commands to use the programs look like.

I will be running the pipeline over the next few days and I'll send you the results so that you can get an idea of resource usage.

## Checking job status

You can check the status of your jobs with `qstat -u <your_username>`

You can get further details for a specific job using `checkjob -v <job_id>`

# Doing a denovo Trinity assembly
Below I will lay out the steps of doing a denovo assembly.

## pre-processing
Before assembly can be undertaken the sequencing data must first go through some QC. Namely, adapter sequences must be checked for and removed and low quality base calls at either end of the sequences should be removed.

Many people like to use [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). I tend to use [fastp](https://github.com/OpenGene/fastp) for its simplicity and speed.

fastp can be [installed](https://anaconda.org/bioconda/fastp) with conda.

You can concatenate all of your sequence files together (i.e. concatenate all R1 files, and concatenate all R2 files) so that you end up with a single R1 and R2 pair. Make sure that you concatenate the files together in the same order.

## Error correction of kmers
Next, kmer correction should be conducted. For RNA-seq data, which has very uneven read distributions, we use [rCorrector](https://github.com/mourisl/Rcorrector).

Unfixable reads need to be discarded manually (see section 4 [here](https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html). I have modified their python script so that it runs with python3. It is here: [remove_unfixable.py](bin/remove_unfixable.py).

## Removal of rRNA reads
This can be acheived by mapping the reads to an rRNA database such as the SILVA database. Follow the advise of the Harvard guide and download the SSUParc and LSUParc fasta files and concatenate them. The files can be downloaded from [here](https://www.arb-silva.de/no_cache/download/archive/current/Exports/). Have a look at the README.txt file in this directory to get an idea of what each of the files contain. I would recommend you use the [SILVA_138.1_LSUParc_tax_silva_trunc.fasta.gz](https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_LSUParc_tax_silva_trunc.fasta.gz) and	[SILVA_138.1_SSUParc_tax_silva_trunc.fasta.gz](https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSUParc_tax_silva_trunc.fasta.gz) files. You can use the `wget` command to download them, and `cat` to concatenate them.

Its best to be conservative and only keep those reads in which neither of the pair mapped to the SILVA db reads.

Building the bowtie2 index or this database takes some time. It may be best if I provide you with the index. Let's discuss.

## Trinity
At this point you're ready to run Trinity. The biggest difficulty with running Trinity will be to get the resource requirements right.

I'll let you know when I've got a better idea of resource use.
