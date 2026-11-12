# !/usr/bin/env python

import sys, os
from Bio.Seq import Seq

variant_file_path = sys.argv[1]
rc_variant_file_path = sys.argv[2]

with open(variant_file_path, "r") as variant_file, open(rc_variant_file_path, "w") as rc_variant_file:

    rc_variant_file.write(variant_file.readline())

    for line in variant_file:

        rc_variant_file.write(line)

        split_line = line.split("\t")

        spec_kmer = split_line[1]
        alt_kmer = split_line[2]

        rc_spec_kmer = str(Seq(spec_kmer).reverse_complement())
        rc_alt_kmer = str(Seq(alt_kmer).reverse_complement())

        split_line[1] = rc_spec_kmer
        split_line[2] = rc_alt_kmer

        rc_variant_file.write("\t".join(split_line))