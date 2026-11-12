import sys

matrix_path = sys.argv[1]
output_path = sys.argv[2] #output variable position file with strain annotated and only containing variable positions unique to one strain

strain_list = []

with open(matrix_path, "r") as matrix_file, open(output_path, "w") as output_file:

    raw_strain_list = matrix_file.readline().strip().split("\t")

    for strain in raw_strain_list:

        cleaned_strain = strain.split("_")[0]
        strain_list.append(cleaned_strain)

    for line in matrix_file:

        split_line = line.strip().split("\t")

        position = split_line[0]

        matrix_line = split_line[1:]

        strains_and_matrix_dict = dict(zip(strain_list, matrix_line))

        strains_and_matrix_dict_values = list(strains_and_matrix_dict.values())

        if strains_and_matrix_dict_values.count("0") == 1:

            for strain in strains_and_matrix_dict:

                if strains_and_matrix_dict[strain] == "0":

                    output_file.write(position + "\t" + strain + "\n")

        elif strains_and_matrix_dict_values.count("1") == 1:

            for strain in strains_and_matrix_dict:

                if strains_and_matrix_dict[strain] == "1":

                    output_file.write(position + "\t" + strain + "\n")

