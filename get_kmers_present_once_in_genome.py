# !/usr/bin/env python

import sys, math, os
from Bio.Seq import Seq

def get_chr_and_pos(chromosome_position_string):

    chr_pos_split = chromosome_position_string.split(":")
    
    return chr_pos_split[0], int(chr_pos_split[1])

fasta_file_path = sys.argv[1]
variable_positions_path = sys.argv[2] #uses the one from get_kmers_dir that has the alt allele already
kmer_freq_path = sys.argv[3]

# make dictionary of variable positions
variable_positions_dict = {}
alt_allele_dict = {}
strain_set = set()

with open(variable_positions_path, "r") as variable_positions_file: 

    for line in variable_positions_file:

        split_line = line.strip().split("\t")
        var_position = split_line[0]
        ref_allele = split_line[1]
        alt_allele = split_line[2]
        strain = split_line[3]
        
        variable_positions_dict[var_position] = strain
        alt_allele_dict[var_position] = alt_allele
        strain_set.add(strain)

print("got reference and alt allele for all variable positions")

strain_list = list(strain_set)

with open(fasta_file_path, "r") as fasta_file:

    fasta_dict = {}
    chromosome = ""
    sequence = ""

    for line in fasta_file:

        #header line
        if line.strip()[0] == ">":

            #add completed chromosome to dict
            if len(chromosome) > 0:

                fasta_dict[chromosome] = sequence

            chromosome = line[1:].split()[0]
            sequence = ""

        else:

            sequence += line.strip()

    #add final chromosome to dict
    fasta_dict[chromosome] = sequence

#print(fasta_dict.keys())
print("got chromosomes")

with open(kmer_freq_path, "w") as kmer_freq_file:

    #get kmer frequencies by length (5 --> 30 bp)
    #for kmer_length in range(5,31):
    kmer_length = 25

    print("looking at kmers of length " + str(kmer_length))

    lower_bound_diff = math.floor(kmer_length/2)

    for variable_position in variable_positions_dict:

        #print("variable position: " + str(variable_position))

        chromosome_name, position = get_chr_and_pos(variable_position)
        alt_allele = alt_allele_dict[variable_position]
        strain = variable_positions_dict[variable_position]

        #print("alt allele: " + alt_allele)
        #print("strain: " + strain)

        lower_bound = position - int(lower_bound_diff) - 1
        upper_bound = lower_bound + kmer_length

        #print("lower_bound: " + str(lower_bound))
        #print("upper_bound: " + str(upper_bound))

        chromosome_sequence = fasta_dict[chromosome_name]

        #print("chromosome name: " + chromosome_name)

        spec_kmer = chromosome_sequence[lower_bound:upper_bound]
        alt_kmer = spec_kmer[:lower_bound_diff] + alt_allele + spec_kmer[lower_bound_diff+1:]

        spec_kmer_rc = str(Seq(spec_kmer).reverse_complement()) 
        alt_kmer_rc = str(Seq(alt_kmer).reverse_complement())

        #print("ref kmer: " + spec_kmer)
        #print("alt kmer: " + alt_kmer)

        spec_freq = 0
        alt_freq = 0
        rc_freq = 0

        #in NF54_3D7, the ref kmer should be in 3D7, while the alt kmer shouldn't be
        if strain == "NF54":

            for chromosome in fasta_dict:

                spec_freq += fasta_dict[chromosome].count(spec_kmer)
                alt_freq += fasta_dict[chromosome].count(alt_kmer)
                rc_freq += fasta_dict[chromosome].count(spec_kmer_rc)
                rc_freq += fasta_dict[chromosome].count(alt_kmer_rc)

            if spec_freq == 1 and alt_freq == 0 and rc_freq == 0:

                kmer_freq_file.write(variable_position + "\t" + spec_kmer + "\t" + alt_kmer + "\t" + strain + "\n")

        #in non NF54_3D7 strains, the alt kmer should be in 3D7, while the ref kmer shouldn't be
        else:

            for chromosome in fasta_dict:

                spec_freq += fasta_dict[chromosome].count(alt_kmer)
                alt_freq += fasta_dict[chromosome].count(spec_kmer)
                rc_freq += fasta_dict[chromosome].count(spec_kmer_rc)
                rc_freq += fasta_dict[chromosome].count(alt_kmer_rc)

            if alt_freq == 1 and spec_freq == 0 and rc_freq == 0:

                kmer_freq_file.write(variable_position + "\t" + spec_kmer + "\t" + alt_kmer + "\t" + strain + "\n")
