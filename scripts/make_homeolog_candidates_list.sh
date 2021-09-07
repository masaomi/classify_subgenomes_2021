#!/bin/sh
# encoding: utf-8

################################################################
# Used software (version)
# NCBI blast (ver. 2.2.31+)
# Ruby (ver. 2.2.3p173)
# BioRuby (RubyGem, ver. 1.5.1)
################################################################

FASTA=references/2Es_HS_3.2.1_transcripts.fa
NUM_THREADS=8

cat <<EOF> estimate_homeologs.rb
#!/usr/bin/env ruby
# encoding: utf-8

ALIGN_LEN = 300

require 'bio'
gene2len = Hash[*Bio::FlatFile.open("$FASTA").to_a.map{|e| [e.definition.split.first, e.seq.length]}.flatten]
all_genes = gene2len.keys
genes = all_genes.length

def blast_besthits
  query = nil
  targets = []
  while line=DATA.gets
    new_query, target, xxx, alignment_length, *others = line.chomp.split
    if query and query != new_query and !targets.empty?
      yield [query, targets.uniq.join(",")]
      query = nil
      targets = []
    end
    if alignment_length.to_i > ALIGN_LEN
      unless query
        query = new_query
      end
      targets << target
    end
  end
  if query and !targets.empty?
    yield [query, targets.uniq.join(",")]
  end
end

uniqs = {}
obhs = {}
blast_besthits do |query, targets|
  targets_ = targets.split(',')
  targets_.delete(query)
  obhs[query] = targets_
  uniqs[query] = targets_.clone
end
duplicates = obhs.keys.length
puts "# Total genes: #{genes}"
puts "# Duplicates (evalue>1.0e-10, alignment length>#{ALIGN_LEN}): #{duplicates}"

uniqs.keys.each do |query|
  if uniqs[query]
    uniqs[query].each do |target|
      uniqs.delete(target)
    end
  end
end
uniq_genes = uniqs.keys.length + (genes - duplicates)
puts "# Unique genes (count duplicates as one): #{uniq_genes}"

dobhs = {}
obhs.each do |query, targets|
  query_length = gene2len[query]
  real_target = nil
  targets.each do |target|
    scaff1 = if query =~ /(Super-Scaffold_\d+)/
               \$1
             elsif query =~ /(scaffold.+_obj)/
               \$1
             end
    scaff2 = if target =~ /(Super-Scaffold_\d+)/
               \$1
             elsif target =~ /(scaffold.+_obj)/
               \$1
             end
    target_length = gene2len[target]
    diff_length = (query_length - target_length).abs.to_f
    ave_length = (query_length + target_length)/2.0
    if scaff1 != scaff2 and diff_length/ave_length < 0.5
      real_target = target
      break
    end
  end
  if real_target and (real_target != query)
    dobhs[query] = real_target
  end
end

rbhs = {}
dobhs.each do |query, target|
  if dobhs[target] == query
    rbhs[query] = target
  end
end
puts "# Reciprocal blast best hits in different scaffolds (Homeologs): #{rbhs.keys.length}"
puts "# Homeologs/genes: #{"%.2f" % (rbhs.keys.length.to_f/genes).to_s}"

open("homeolog_candidates_list.dat", "w") do |out|
  rbhs.each do |g1, g2|
    out.puts [g1, g2].join("\t")
  end
end
puts
puts "# homeolog_candidates_list.dat generated"

__END__
EOF

echo
echo "# estimate_homeologs.rb generated"
makeblastdb -dbtype nucl -hash_index -in $FASTA
blastn -num_threads $NUM_THREADS -query $FASTA -db $FASTA -evalue 1.0e-10 -outfmt 6 >> estimate_homeologs.rb
echo "# `which blastn | sed -e "s/[\r\n]\+//g"` executed"
echo
echo "# ruby estimate_homeologs.rb"
ruby estimate_homeologs.rb

