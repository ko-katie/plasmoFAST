# !/usr/bin/env python

import sys

original_variable_position_path = sys.argv[1]
filtered_variable_position_path = sys.argv[2]

with open(original_variable_position_path, "r") as original_variable_position_file, open(filtered_variable_position_path, "w") as filtered_variable_position_file:

    first_line = original_variable_position_file.readline()
    filtered_variable_position_file.write(first_line)

    first_line_split = first_line.strip().split("\t")
    prev_pos = int(first_line_split[0].split(":")[-1])
    prev_chrom = first_line_split[0].split(":")[0]

    for line in original_variable_position_file:

        split_line = line.strip().split("\t")

        variable_position_string = split_line[0]

        curr_pos = int(variable_position_string.split(":")[1])
        curr_chrom = variable_position_string.split(":")[0]

        if prev_chrom != curr_chrom or curr_pos - prev_pos >= 35:

            filtered_variable_position_file.write(line)

            prev_pos = curr_pos
            prev_chrom = curr_chrom