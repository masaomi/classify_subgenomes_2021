# Classifying subgenomes for a polyploid assembled genome

This is the script to classify subgenomes for a allo-tetraploid assembled genome based on the alignment results of either one parental sequenced reads.

Requirements
* references/
   * 2Es_HS_3.2.1_genome.fa: assembled result (fasta)
   * 2Es_HS_3.2.1_transcripts.fa: predicted transcripts (fasta)
   * 2Es_HS_3.2.1_maker.gff: predicted annotation
* Eind.sort.bam: alignment result of one parental diploid to the reference

Run
```
$ sh classify_subgenomes.sh
```

