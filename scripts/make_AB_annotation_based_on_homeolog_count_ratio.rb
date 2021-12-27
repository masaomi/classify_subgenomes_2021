#!/usr/bin/env ruby
# encoding: utf-8
# Version = '20211227-140831'

HOMEOLOG_RATIO_THREAHOLD = 0.5
unless genome_fa=ARGV[0] and genes_gff=ARGV[1] and scaff_ab_homeolog_count_ratio_dat=ARGV[2] and homeolog_list_dat=ARGV[3]
  puts <<-eos
  usage:
    #{File.basename(__FILE__)} [genome.fa] [genes.gff] [scaff_ab_homeolog_count_ratio.dat] [homeolog_list.dat]
  option:
    --homeolog_ratio_threshold: the ratio of A-genes in the scaffold (default: #{HOMEOLOG_RATIO_THRESHOLD})
    --out_dir: default .
  eos
  exit
end

homeolog_ratio_threshold = if idx = ARGV.index("--homeolog_ratio_threshold")
                             ARGV[idx+1].to_f
                           else
                             HOMEOLOG_RATIO_THREAHOLD
                           end
out_dir = if idx = ARGV.index("--out_dir")
            require "fileutils"
            FileUtils.mkdir_p ARGV[idx+1]
            ARGV[idx+1]
          else
            "."
          end

require 'csv'
sid2ab = {}
asid = 0
bsid = 0
nsid = 0
CSV.foreach(scaff_ab_homeolog_count_ratio_dat, :headers=>true, :col_sep=>"\t") do |row|
  sid = row["scaffold"]
  a_homeolog_ratio = row["A_homeolog_ratio"].to_f
  b_homeolog_ratio = 1.0 - a_homeolog_ratio
  if a_homeolog_ratio >= homeolog_ratio_threshold
    sid2ab[sid] = "sa#{"%04d" % (asid+=1)}"
  elsif b_homeolog_ratio >= homeolog_ratio_threshold
    sid2ab[sid] = "sb#{"%04d" % (bsid+=1)}"
  else
    sid2ab[sid] = "sn#{"%04d" % (nsid+=1)}"
  end
end
#p sid2ab
#exit

gff_lines = {}
gid2ab = {}
agid = 0
bgid = 0
ngid = 0
File.readlines(genes_gff).each do |line|
  if line =~ /^#/
    gff_lines[:header] ||= []
    gff_lines[:header] << line.chomp
  else
  # Super-Scaffold_5  maker gene  497241  499468  . + . ID=maker-Super-Scaffold_5-augustus-gene-0.35;Name=maker-Super-Scaffold_5-augustus-gene-0.35
  # Super-Scaffold_5  maker mRNA  497241  499468  . + . ID=maker-Super-Scaffold_5-augustus-gene-0.35-mRNA-1;Parent=maker-Super-Scaffold_5-augustus-gene-0.35;Name=maker-Super-Scaffold_5-augustus-gene-0.35-mRNA-1;_AED=0.11;_eAED=0.11;_QI=0|0|0|0.75|0.66|0.75|4|0|309
  # Super-Scaffold_5  maker exon  497241  497489  . + . ID=maker-Super-Scaffold_5-augustus-gene-0.35-mRNA-1:exon:15696;Parent=maker-Super-Scaffold_5-augustus-gene-0.35-mRNA-1
  # Super-Scaffold_5  maker CDS 497241  497489  . + 0 ID=maker-Super-Scaffold_5-augustus-gene-0.35-mRNA-1:cds;Parent=maker-Super-Scaffold_5-augustus-gene-0.35-mRNA-1
  # Super-Scaffold_5  maker five_prime_UTR  11866356  11866444  . + . ID=maker-Super-Scaffold_5-augustus-gene-11.101-mRNA-1:five_prime_utr;Parent=maker-Super-Scaffold_5-augustus-gene-11.101-mRNA-1
    sid_, maker, type, *others = line.split
    sid = sid_.gsub('|', '_').gsub('_obj', '')
    if type == "gene" and line =~ /ID=(.*?);/
      gid = $1
      case sid2ab[sid]
      when /sa/
        #gid2ab[gid] = "ga#{"%05d" % (agid+=1)}"
        gid2ab[gid] = "ga"
      when /sb/
        #gid2ab[gid] = "gb#{"%05d" % (bgid+=1)}"
        gid2ab[gid] = "gb"
      else
        sid2ab[sid] ||= "sn#{"%04d" % (nsid+=1)}"
        #gid2ab[gid] = "gn#{"%05d" % (ngid+=1)}"
        gid2ab[gid] = "gn"
      end
    end
    gff_lines[sid] ||= []
    gff_lines[sid] << line.chomp
  end
end
#p gid2ab
#p gff_lines
#puts gff_lines["sa0001"].join("\n")
#exit

# new gids
scaff_pos_gids = {}
sid2ab.sort_by{|osid, nsid| nsid}.each do |osid, nsid|
  gff_lines[osid].each do |line|
    sid, maker, type, sp, ep, dot1, strand, dot2, *others = line.split
    raise unless others.first =~ /ID=/
    if type == "gene"
      scaff_pos_gids[nsid] ||= []
      gid = others.join.match(/ID=(.*?);/)[1].split(/-mRNA-1/).first
      scaff_pos_gids[nsid] << [sp.to_i, gid]
    end
  end
end
#p scaff_pos_gids
#exit
ogid2ngid = {}
sid2ab.sort_by{|osid, nsid| nsid}.each do |osid, nsid|
  scaff_pos_gids[nsid].sort_by{|sp, gid| sp}.each do |sp, gid|
    gtype = gid2ab[gid]
    new_gid = case gtype
              when "ga"
                "ga#{"%05d" % (agid+=1)}"
              when "gb"
                "gb#{"%05d" % (bgid+=1)}"
              when "gn"
                "gn#{"%05d" % (ngid+=1)}"
              else
                raise
              end
    ogid2ngid[gid] = new_gid
  end
end
#p ogid2ngid
#exit

exon_no = -1
new_genes_gff = File.join(out_dir, File.basename(genes_gff).gsub(/.gff/, '_new.gff'))
open(new_genes_gff, "w") do |out|
  out.puts gff_lines[:header].join("\n")
  sid2ab.sort_by{|osid, nsid| nsid}.each do |osid, nsid|
    begin
      gff_lines[osid].each do |line|
        begin
          sid, maker, type, sp, ep, dot1, strand, dot2, *others = line.split
          raise unless others.first =~ /ID=/
          new_sid = nsid
          new_row = [nsid, "maker", type, sp, ep, dot1, strand, dot2]
          gid = others.join.match(/ID=(.*?);/)[1].split(/-mRNA-1/).first.split('.').first
          new_gid = ogid2ngid[gid]
          unless new_gid
            p gid
            raise
          end
          case type
          when "gene"
            new_row << "ID=#{new_gid};Name=#{new_gid}" 
          when "mRNA"
            aed = others.join.match(/(_AED.+$)/)[1]
            new_row << "ID=#{new_gid}.t1;Parent=#{new_gid};Name=#{new_gid}.t1;#{aed}" 
          when "exon"
            new_row << "ID=#{new_gid}.t1:exon:#{exon_no+=1};Parent=#{new_gid}.t1"
          when "CDS"
            new_row << "ID=#{new_gid}.t1:cds;Parent=#{new_gid}.t1"
          when "three_prime_UTR"
            new_row << "ID=#{new_gid}.t1:three_prime_utr;Parent=#{new_gid}.t1"
          when "five_prime_UTR"
            new_row << "ID=#{new_gid}.t1:five_prime_utr;Parent=#{new_gid}.t1"
          else
            raise
          end
          out.puts new_row.join("\t")
        rescue
          p line
          exit
        end
      end
    rescue
      p osid
      exit
    end
  end
end
puts "# #{new_genes_gff} generated"

# genome_new.fa
require 'bio'
new_genome_fa = File.join(out_dir, File.basename(genome_fa).gsub(/.fa/, '_new.fa'))
sid2seq = {}
Bio::FlatFile.open(genome_fa).each do |e|
  sid_ = e.definition
  sid = sid_.gsub('|', '_').gsub('_obj', '')
  sid2ab[sid] ||= "sn#{"%04d" % (nsid+=1)}"
  sid2seq[sid] = e.seq
end
open(new_genome_fa, "w") do |out|
  sid2ab.sort_by{|osid, nsid| nsid}.each do |osid, nsid|
    out.puts ">#{nsid}"
    out.puts sid2seq[osid].scan(/.{100}|.+\Z/).join("\n")
  end
end
puts "# #{new_genome_fa} generated"

require 'fileutils'
unless File.dirname(homeolog_list_dat) == out_dir
  FileUtils.cp homeolog_list_dat, out_dir
  puts "# #{homeolog_list_dat} copied to #{out_dir}"
end
unless File.dirname(scaff_ab_homeolog_count_ratio_dat) == out_dir
  FileUtils.cp scaff_ab_homeolog_count_ratio_dat, out_dir
  puts "# #{scaff_ab_homeolog_count_ratio_dat} copied to #{out_dir}"
end
