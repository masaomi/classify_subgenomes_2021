#!/usr/bin/env ruby
# encoding: utf-8
# Version = '20181128-083523'

unless gff=ARGV[0] and sorted_bam=ARGV[1] and out_dir=ARGV[2]
  puts "usage:"
  puts " ruby #{File.basename(__FILE__)} [gff] [sorted_bam] [out_dir] > batch.sh"
  exit
end

require 'fileutils'

FileUtils.mkdir_p out_dir

genes = {}
File.readlines(gff).each do |line|
  unless line =~ /#/
    scaff, maker, type, sp, ep, *others = line.split
    if type == "gene" and line =~ /ID=(.+?);/
      gid = $1
      genes[gid] = [scaff, sp.to_i, ep.to_i]
    end
  end
end
genes.each do |gid, ginfo|
  chr, sp, ep = ginfo
  command = "samtools depth -aa -r #{chr.gsub('|', '\|')}:#{sp}-#{ep} #{sorted_bam} > #{out_dir}/#{gid.gsub('|', '_')}.depth"
  puts command
end
