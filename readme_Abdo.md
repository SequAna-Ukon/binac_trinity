# down-stream Trinity de-novo assembly on BINAC

## Transcriptome Assembly Quality Assessment


- Examine the RNA-Seq read representation of the assembly. Ideally, at least ~80% of your input RNA-Seq reads are represented by your transcriptome assembly. The remaining unassembled reads likely corresponds to lowly expressed transcripts with insufficient coverage to enable assembly, or are low quality or aberrant reads.

- Examine the representation of full-length reconstructed protein-coding genes, by searching the assembled transcripts against a database of known protein sequences.

- Use BUSCO to explore completeness according to conserved ortholog content.

- Compute the E90N50 transcript contig length - the contig N50 value based on the set of transcripts representing 90% of the expression data.

- Compute DETONATE scores. DETONATE provides a rigorous computational assessment of the quality of a transcriptome assembly, and is useful if you want to run several assemblies using different parameter settings or using altogether different tools. That assembly with the highest DETONATE score is considered the best one.

- Try using TransRate. TransRate generates a number of useful statistics for evaluating your transcriptome assembly. Read about TransRate here: http://genome.cshlp.org/content/26/8/1134. Note that certain statistics may be biased against the large numbers of transcripts that are very lowly expressed. Consider generating TransRate statistics for your transcriptome before and after applying a minimum expression-based filter.

- Explore rnaQUAST a quality assessment tool for de novo transcriptome assemblies.
