# !/usr/bin/env python3

import sys, time
import pandas as pd
import matplotlib.pyplot as plt

#determines category of variable position based on specific and nonspecific coverage
def get_category(spec_cov, nonspec_cov, variable_position_name, strain):

    total_coverage = spec_cov + nonspec_cov

    #if fewer than 30 occurences of kmers for position are found, return low coverage
    if total_coverage < 30:

        return 3

    #more than 90% of kmers are specific, return specific category
    elif spec_cov / (total_coverage) >= 0.90:

        return 0

    #less than 10% of kmers are specific, return nonspecific category
    elif spec_cov / (total_coverage) <= 0.10:

        return 1

    #between 10% to 90% of kmers are specific, return mixed category
    else:

        return 2

start_time = time.time()

#path to file containg variable positions and their kmers
kmer_list_path = sys.argv[1]
#path to file containing kmc output of all kmers found in sample and their counts
kmc_output_path = sys.argv[2]
#path to write output file containing number of spec, nonspec, mixed, and low coverage positions per strain
output_file_path = sys.argv[3]
#path to barplot of positions
barplot_path = sys.argv[4]

#store variable positions as keys, and list as value: [spec_cov, nonspec_cov, strain]
variable_position_dict = {}
#store specific kmers as keys and their variable position as value
spec_kmer_dict = {}
#store nonspecific kmers as keys and their variable position as value
nonspec_kmer_dict = {}

#open kmer file and read header line
kmer_list_file = open(kmer_list_path, 'r')
kmer_list_file.readline()

#read through kmer_list file
for line in kmer_list_file:

    #split tab-delimited line
    split_line = line.strip().split("\t")
    
    chromosome = split_line[0]
    position = split_line[1]
    strain_specific_kmer = split_line[2]
    nonspecific_kmer = split_line[3]
    strain = split_line[4]

    #build variable position name as "chromosome:position"
    variable_position_name = chromosome + ":" + position

    #add variable position to dictinoary, initiliaze value with [freq of specific kmer, freq of nonspecific kmer, strain]
    variable_position_dict[variable_position_name] = [0, 0, strain]
    #add specific and nonspecific kmers to their respective dictionaries
    spec_kmer_dict[strain_specific_kmer] = variable_position_name
    nonspec_kmer_dict[nonspecific_kmer] = variable_position_name

kmer_list_file.close()

#open kmc output file
kmc_output_file = open(kmc_output_path, 'r')

#track current kmer and its count for reading through kmc output file
curr_kmer = ""
curr_count = 0

#read through kmc output file
for line in kmc_output_file:

    #split tab-delimited line
    split_line = line.strip().split("\t")

    #get current kmer and its count
    curr_kmer = split_line[0]
    curr_count = int(split_line[1])
            
    #check if speciic kmer
    if curr_kmer in spec_kmer_dict:

        #add to specific kmer dict and increase specific kmer count for variable position by curr_count
        variable_position_name = spec_kmer_dict[curr_kmer]
        variable_position_dict[variable_position_name][0] += curr_count

    #check if nonspecific kmer
    elif curr_kmer in nonspec_kmer_dict:

        #add to nonspecific kmer dict and increase nonspecific kmer count for variable position by curr_count
        variable_position_name = nonspec_kmer_dict[curr_kmer]
        variable_position_dict[variable_position_name][1] += curr_count

kmc_output_file.close()

#initialize strain_output_dict to keep track of number of spec, nonspec, mix, and low coverage positions per strain
strain_list = ["NF54_3D7", "7G8", "W2_Dd2", "D10", "HB3", "GB4"]
strain_output_dict = dict.fromkeys(strain_list)

#initialize all types of positions for each strain with 0
for strain in strain_output_dict:

    strain_output_dict[strain] = [0,0,0,0]

#go through all variable positions and update strain_output_dict based on coverage and allele frequency
for variable_position_name in variable_position_dict:

    #get variable position information from variable_position_dict
    variable_position = variable_position_dict[variable_position_name]

    spec_cov = variable_position[0]
    nonspec_cov = variable_position[1]
    strain = variable_position[2]
    
    #update strain_output_dict for the strain of the variable position based on category determined by get_category function
    strain_output_dict[strain][get_category(spec_cov, nonspec_cov, variable_position_name, strain)] += 1

#open output file
output_file = open(output_file_path, "w")

#write header
output_file.write("Strain\tSpec Positions\tNonspec Positions\tMixed Positions\tLow Coverage Positions\n")

#write output for each strain to file
for strain in strain_output_dict:

    print(str(strain_output_dict[strain]))

    output_file.write(strain + "\t" + "\t".join(map(str, strain_output_dict[strain])) + "\n")

output_file.close()

#print run time
end_time = time.time()
print("time ran: " + str(end_time - start_time))

### PLOT RESULTS ###
data=[strain_output_dict[strain] for strain in strain_list]

#set strain names as the index
df = pd.DataFrame(data, index=strain_list, columns=['Specific', 'Nonspecific', 'Mixed', 'Low Coverage Positions'])

#df = df.drop('Low Coverage Positions', axis=1)

#plot stacked bar chart
ax = df.plot(
    kind="bar",
    stacked=True,
    figsize=(8,6),
    color=["#31a087",'#e64b35','#3c5488', "gray"]
)

#labels and formatting
ax.set_xlabel("Strain")
ax.set_ylabel("Number of Positions")
ax.set_title("Position Types per Strain")
plt.xticks(rotation=45)
plt.tight_layout()

#save and show
plt.savefig(barplot_path, dpi=300)