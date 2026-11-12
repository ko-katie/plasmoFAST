# !/usr/bin/env python3

import sys, os
from Bio import SeqIO
from Bio.Seq import Seq

kmer_list_path = sys.argv[1]
filtered_kmer_list_path = sys.argv[2]
human_genome_path = sys.argv[3]

kmer_dict = {}
variable_position_dict = {}

with open(kmer_list_path, "r") as kmer_list_file:

    kmer_list_file.readline()

    for line in kmer_list_file:

        split_line = line.strip().split("\t")
        variable_position = split_line[0]
        spec_kmer = split_line[1]
        alt_kmer = split_line[2]
        strain = split_line[3]

        rc_spec_kmer = str(Seq(spec_kmer).reverse_complement())
        rc_alt_kmer = str(Seq(alt_kmer).reverse_complement())

        kmer_dict[spec_kmer] = variable_position
        kmer_dict[alt_kmer] = variable_position
        kmer_dict[rc_spec_kmer] = variable_position
        kmer_dict[rc_spec_kmer] = variable_position

        variable_position_dict[variable_position] = 0

print(kmer_dict)
with open(human_genome_path) as human_genome_handler:

    for record in SeqIO.parse(human_genome_handler, "fasta"):

        if record.id[:2] == "NC":

            print("record.id: " + record.id)
            sequence = record.seq.upper()

            n = 0

            while n <= len(sequence) - 25:

                lower_bound = n
                upper_bound = n+25

                string_slice = sequence[lower_bound:upper_bound]
                #print("string slice: " + string_slice)

                if string_slice in kmer_dict:

                    print("string slice found: " + string_slice, flush=True)
                    variable_position = kmer_dict[string_slice]

                    variable_position_dict[variable_position] += 1
            
                n += 1

with open(filtered_kmer_list_path, "w") as filtered_kmer_list_file, open(kmer_list_path, "r") as kmer_list_file:

    kmer_list_file.readline()

    for line in kmer_list_file:

        split_line = line.split("\t")
        variable_position = split_line[0]

        if variable_position_dict[variable_position] == 0:

            filtered_kmer_list_file.write(line)


