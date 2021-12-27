#!/usr/bin/env ruby
# encoding: utf-8
# Version = '20211227-115754'

require 'csv'

MIN_COV = 20.0 # minimal coverage to be A genome
MAX_RATIO = 0.2

#unless min_cov=ARGV[0] and max_ratio=ARGV[1]
unless genes_gff=ARGV[0] and gene_depth_length_dat=ARGV[1] and homeolog_candidates_list_dat=ARGV[2]
  puts "usage:"
  puts " ruby #{File.basename(__FILE__)} [genes.gff] [gene_depth_length.dat] [homeolog_candidates_list.dat]"
  puts
  puts "option:"
  puts " --batch: on (default: off)"
  puts " --min_cov: min_cov (default: #{MIN_COV})"
  puts " --max_ratio: max_ratio (default: #{MAX_RATIO})"
  exit
end

batch_mode = ARGV.index("--batch")
min_cov = if idx = ARGV.index("--min_cov")
            ARGV[idx+1].to_f
          else
            MIN_COV
          end
max_ratio = if idx = ARGV.index("--max_ratio")
              ARGV[idx+1].to_f
            else
              MAX_RATIO
            end

# load gene depth
gene2depth = {}
gene2length = {}
CSV.foreach(gene_depth_length_dat, :headers=>true, :col_sep=>"\t") do |row|
  gid = row["GID"]
  dep = row["Depth_Ave"].to_f
  len = row["Length"].to_i
  gene2depth[gid] = dep
  gene2length[gid] = len
end

# load gff
rid2len = {}
rid2scaff = {}
File.readlines(genes_gff).each do |line|
  scaffid, software, type, sp, ep, *others = line.chomp.split
  if type == "mRNA" and line =~ /ID=(.+);Parent/
    rid = $1
    rid2len[rid] = (ep.to_i - sp.to_i).abs + 1
    rid2scaff[rid] = scaffid
  end
end

# A B identifiy
homeolog_count = 0
out = nil
unless batch_mode
  out = open("homeolog_list.dat", "w") 
  out.puts ["A_homeolog", "B_homeolog", "A_coverage", "B_coverage", "coverage_ave", "A_coverage_ratio", "A_length", "B_length", "length_ave", "length_ratio(A/B)", "A_scaffold", "B_scaffold"].join("\t")
end
already_out = {}
length_ratios = []
File.readlines(homeolog_candidates_list_dat).each do |line|
  r1, r2 = line.chomp.split
  s1 = rid2scaff[r1]
  s2 = rid2scaff[r2]
  g1 = r1.split('-mRNA').first.gsub('|', '_').gsub(/.t1$/, '')
  g2 = r2.split('-mRNA').first.gsub('|', '_').gsub(/.t1$/, '')
  if g1d = gene2depth[g1] and g2d = gene2depth[g2] and
     g1l = gene2length[g1] and g2l = gene2length[g2]
    if !already_out[g1] and !already_out[g2]
      min, max = [g1d, g2d].sort 
      sum = min + max
      if max > min_cov
        homeolog_ratio = min / sum
        if homeolog_ratio < max_ratio
          homeolog_count += 1
          unless batch_mode
            length_ave = (g1l+g2l)/2.0
            coverage_ave = (g1d*g1l+g2d*g2l)/(g1l+g2l)
            if g1d > g2d
              length_ratio = g1l/g2l.to_f
              length_ratios << length_ratio
              a_coverage_ratio = g1d/(g1d+g2d)
              out.puts [g1, g2, g1d, g2d, coverage_ave, a_coverage_ratio, g1l, g2l, length_ave, length_ratio, s1, s2].join("\t")
            else
              length_ratio = g2l/g1l.to_f
              length_ratios << length_ratio
              a_coverage_ratio = g2d/(g1d+g2d)
              out.puts [g2, g1, g2d, g1d, coverage_ave, a_coverage_ratio, g2l, g1l, length_ave, length_ratio, s2, s1].join("\t")
            end
          end
        end
      end
      already_out[g1] = true
      already_out[g2] = true
    end
  end
end
unless batch_mode
  out.close
end

unless batch_mode
  puts "# homeolog_list.dat generated"
  puts
  puts "# MIN_COV: #{min_cov}"
  puts "# MAX_RATIO: #{max_ratio}"
  puts "# detected homeolog count: #{homeolog_count}*2 = #{homeolog_count*2}"
  puts "# max gene length ratio: #{length_ratios.max}"
  puts "# min gene length ratio: #{length_ratios.min}"
else
  puts [min_cov, max_ratio, homeolog_count*2].join("\t")
end
