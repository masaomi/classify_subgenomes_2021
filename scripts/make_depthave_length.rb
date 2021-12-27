#!/usr/bin/env ruby
# encoding: utf-8
# Version = '20211227-112302'

unless depth_dir=ARGV[0] and transcripts_fa=ARGV[1] and gff=ARGV[2]
  puts "usage:"
  puts " ruby #{File.basename(__FILE__)} [depth_dir] [transcripts_fa] [gff] > gene_depth_length.dat"
  exit
end

require 'bio'
trans2len = Hash[*Bio::FlatFile.open(transcripts_fa).to_a.map{|e| [e.definition.split.first.split('-mRNA').first.gsub('|', '_'), e.seq.length]}.flatten]

gene2len = {}
File.readlines(gff).each do |line|
  # Super-Scaffold_1  maker gene  230988  234682  . + . ID=maker-Super-Scaffold_1-augustus-gene-0.136;Name=maker-Super-Scaffold_1-augustus-gene-0.136
  scaffid, software, type, sp, ep, *others = line.chomp.split
  if type == "gene" and line =~ /ID=(.+);Name/
    gid = $1
    gid = gid.gsub('|', '_')
    gene2len[gid] = (ep.to_i - sp.to_i).abs + 1
  end
end

samtools_len_transcript_len = []
samtools_len_gene_len = []
puts ["GID", "Depth_Ave", "Length", "Gene_Length", "Transcript_Length"].join("\t")
Dir[File.join(depth_dir, "*.depth")].sort.each do |file|
  gid = File.basename(file).gsub(".depth", "")
  count = 0
  ave_depth = 0
  File.readlines(file).each do |line|
    scaff, pos, cov = line.chomp.split
    ave_depth += cov.to_i
    count += 1
  end
  ave_depth /= count.to_f
  gene_length = gene2len[gid]
  tid = gid + ".t1"
  tran_length = (trans2len[gid]||trans2len[tid])
  unless gene_length and tran_length
    raise gid
  end
  samtools_len_transcript_len << count.to_f/tran_length
  samtools_len_gene_len << count.to_f/gene_length
  puts [gid, ave_depth, count, gene_length, tran_length].join("\t")
end
warn "# 1. max samtools_len_vs_transcript_len: #{samtools_len_transcript_len.max}"
warn "# 2. min samtools_len_vs_transcript_len: #{samtools_len_transcript_len.min}"
warn "# 3. max samtools_len_vs_gene_len: #{samtools_len_gene_len.max}"
warn "# 4. min samtools_len_vs_gene_len: #{samtools_len_gene_len.min}"
warn "# if (3)==(4)==1.0, which means samtools_len == gene_len, and Depth_Ave is calculated by gene length including exons and introns"

