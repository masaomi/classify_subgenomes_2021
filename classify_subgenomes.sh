#!/usr/bin/bash
# encoding: utf-8
# Version = '20211227-133746'

# required
# * references/
#    2Es_HS_3.2.1_genome.fa: assembled result (fasta)
#    2Es_HS_3.2.1_transcripts.fa: predicted transcripts (fasta)
#    2Es_HS_3.2.1_maker.gff: predicted annotation
# * Eind.sort.bam: alignment result of one parental diploid to the reference

sh make_homeolog_candidates_list.sh 
# Total genes: 62554
# Duplicates (evalue>1.0e-10, alignment length>300): 57532
# Unique genes (count duplicates as one): 34101
# Reciprocal blast best hits in different scaffolds (Homeologs): 33786
# Homeologs/genes: 0.54
# homeolog_candidates_list.dat generated

samtools index bwa_out_Ei_350_on_2Es_HS_3.2.1_genome/Eind.sort.bam
ruby scripts/make_samtools_depth.rb references/2Es_HS_3.2.1_maker.gff Ei_350.sort.bam samtools_depth_out > samtools_depth.sh

sh samtools_depth.sh

ruby scripts/make_depthave_length.rb samtools_depth_out references/2Es_HS_3.2.1_transcripts.fa references/2Es_HS_3.2.1_maker.gff > gene_depth_length.dat
# 1. max samtools_len_vs_transcript_len: 31.109929078014183
# 2. min samtools_len_vs_transcript_len: 1.0
# 3. max samtools_len_vs_gene_len: 1.0
# 4. min samtools_len_vs_gene_len: 1.0
# if (3)==(4)==1.0, which means samtools_len == gene_len, and Depth_Ave is calculated by gene length including exons and introns

ruby scripts/make_homeolog_list.rb references/2Es_HS_3.2.1_maker.gff gene_depth_length.dat homeolog_candidates_list.dat
# homeolog_list.dat generated

ruby scripts/scaff_ab_homeolog_count_ratio.rb references/2Es_HS_3.2.1_genome.fa references/2Es_HS_3.2.1_maker.gff homeolog_list.dat --a_homeolog_ratio 0.5 2> /dev/null
# scaff_ab_homeolog_count_ratio.dat generated

ruby scripts/make_AB_annotation_based_on_homeolog_count_ratio.rb references/2Es_HS_3.2.1_genome.fa references/2Es_HS_3.2.1_maker.gff scaff_ab_homeolog_count_ratio.dat homeolog_list.dat --homeolog_ratio_threshold 0.6 --out_dir 2Es_classify_out_homeolog_ratio_0.6

