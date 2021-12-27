#!/usr/bin/env ruby
# encoding: utf-8
# Version = '20211227-133053'

require 'csv'

A_HOMEOLOG_RATIO = 0.5

unless genome_fa=ARGV[0] and genes_gff=ARGV[1] and homeolog_list_dat=ARGV[2]
  puts <<-eos
  usage:
    #{File.basename(__FILE__)} [genome.fa] [genes.gff] [homeolog_list.dat]
  note:
    you can skip nil case by 2> /dev/null
  eos
  exit
end

# load scaffold length
require 'bio'
scaff2len = {}
Bio::FlatFile.open(genome_fa).each do |e|
  sid = e.definition.gsub('|', '_').gsub('_obj', '')
  scaff2len[sid] = e.seq.length
end
#p scaff2len
#exit

gid2sid = {}
File.readlines(genes_gff).each do |line|
  unless line =~ /^#/
    sid, maker, type, *others = line.split
    if type == "gene" and line =~ /ID=(.+?);/
      gid = $1
      gid2sid[gid] = sid
    end
  end
end

scaff_a_count = {}
scaff_b_count = {}
CSV.foreach(homeolog_list_dat, :headers=>true, :col_sep=>"\t") do |row|
  a_gene = row["A_homeolog"]
  b_gene = row["B_homeolog"]
  scaff_a_gene = gid2sid[a_gene]
  scaff_b_gene = gid2sid[b_gene]
  scaff_a_count[scaff_a_gene] ||= 0
  scaff_a_count[scaff_a_gene] += 1
  scaff_b_count[scaff_b_gene] ||= 0
  scaff_b_count[scaff_b_gene] += 1
end

all_scaff = scaff_a_count.keys.concat(scaff_b_count.keys)
all_scaff = all_scaff.uniq.sort
#p all_scaff
#exit

sid2ab = {}
open("scaff_ab_homeolog_count_ratio.dat", "w") do |out|
  out.puts ["scaffold", "A_gene_count", "B_gene_count", "A_homeolog_ratio", "scaff_len"].join("\t")
  all_scaff.sort_by{|scaff| scaff2len[scaff]}.reverse.each do |scaff|
    a_count = scaff_a_count[scaff].to_i
    b_count = scaff_b_count[scaff].to_i
    a_homeolog_ratio = a_count.to_f/(a_count+b_count)
    b_homeolog_ratio = 1.0 - a_homeolog_ratio.to_f
    scaff_len = scaff2len[scaff]

    unless scaff_len
      warn scaff
      raise
    end

    if a_homeolog_ratio > A_HOMEOLOG_RATIO
      sid2ab[scaff] = "A"
    else # b_homeolog_ratio >= A_HOMEOLOG_RATIO
      sid2ab[scaff] = "B"
    end

    a_homeolog_ratio = "%.3f" % a_homeolog_ratio
    out.puts [scaff, a_count, b_count, a_homeolog_ratio, scaff_len].join("\t")
  end
end
puts "# scaff_ab_homeolog_count_ratio.dat generated"

# check
a_scaffold = sid2ab.values.count("A")
b_scaffold = sid2ab.values.count("B")
#u_scaffold = sid2ab.values.count(nil)
total = sid2ab.keys.length

# load homeolog list
require 'csv'
count = 0
success = 0
failure = 0
total_homeolog = 0
CSV.foreach(homeolog_list_dat, :headers=>true, :col_sep=>"\t") do |row|
  total_homeolog += 1
  a_homeolog = row["A_homeolog"]
  b_homeolog = row["B_homeolog"]
  a_sid = row["A_scaffold"]
  b_sid = row["B_scaffold"]
  a_subgenome = sid2ab[a_sid]
  b_subgenome = sid2ab[b_sid]
  if a_subgenome and b_subgenome
    unless a_subgenome == b_subgenome
      success += 1
    else
      failure += 1
    end
    count += 1
  else
  warn ["nil case", a_sid, b_sid].join("\t")
  end
end
puts
puts "# A-subgenome homeolog ratio threshold: #{A_HOMEOLOG_RATIO} (assuming)"
puts "# A_scaffold: #{a_scaffold}, B_scaffold: #{b_scaffold}, total scaffold: #{total}"
puts "# Total homeolog pairs: #{total_homeolog}"
puts "# Total comparisons: #{count} (success:#{success}, failure:#{failure})"
puts "# Total comparisons/Total homeolog pairs = #{"%.2f" % (count/total_homeolog.to_f)}"
puts "# Homeolog separation ratio = success/total comparisons = #{"%.2f" % (success.to_f/count)}"
