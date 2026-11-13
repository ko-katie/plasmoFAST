# plasmoFAST
k-mer based tool for detecting Plasmodium falciparum lab strains from sequencing data

Developed using [kmc (v3.2.4)](https://github.com/refresh-bio/KMC), python3 (v3.8.2), pandas, and matplotlib

### File Overview

**run_kmc_new.sh:** Bash script used to run plasmoFAST. First runs kmc on sequencing data and identify frequency of variable positions. Kmc generates a list of all kmers found in sequencing data and their frequency, which is passed to parse_kmc_output.py to determine the frequency of variable positions.

```
bash run_kmc_new.sh /path/to/working_directory /path/to/sample_input_file.txt
```
- **/path/to/working_directory** is the path to the working directory for analysis to be run in (note lack of trailing backslash)
- **/path/to/sample_input_file.txt** is the path to a txt file containing path(s) to fastq file(s) for sample, one path on each line. Sample name is extracted from file name 

**parse_kmc_output.py:** Python script called by run_kmc_new.sh to parse kmc output and determine frequency of variable positions. Script first reads in kmers for variable positions provided in 25mer_rc_list.txt and stores in dictionaries based on whether they are strain specific or non-strain specific. Next, the output of kmc is read in and it determines the frequency of strain specific and non-specific kmers for each variable position. Based on the number of kmers found for each variable position, and the proportion of specific vs non-specific positions, plasmoFAST determines if each variable position is "Low Coverage", "Specific", "Nonspecific", or "Mixed". Pandas is then used to provide a stacked barplot of, for each strain, the proportion of positions found. 

**25mer_rc_list.txt:** Tab-delimited text file containing curated list of 25mers used to differentiate Pf lab strains. For each variable position plasmoFAST uses, there are two lines, for the forward and reverse sequences of each kmer. Each line contains the chromosome, position of variable position, strain specific 25bp kmer, non-specific 25bp kmer, and strain. This file is read in by parse_kmc_output.py to be stored as a dictionary and searched for in user-provided sequencing data.

### Output
plasmoFAST outputs a {sample_name}_strain_output.txt file with the breakdown of "Low Coverage", "Specific", "Nonspecific", or "Mixed" positions found for each laboratory strain, as well as a stacked barplot representing this data.

# K-mer generation
The following section details the process by which k-mers were selected for use with _PlasmoFAST_
Code developed using: [Hisat2 version 2.2.1](https://daehwankimlab.github.io/hisat2/download/), [gatk version 4.2.2.0](https://github.com/broadinstitute/gatk/releases), [samtools version 1.9](https://www.htslib.org/download/), [R version 4.4.2](https://cran.rstudio.com/), [picard version 2.9.4](https://broadinstitute.github.io/picard/), [vcftools version 0.1.15](https://sourceforge.net/projects/vcftools/), [slurm version 23.11.6](https://slurm.schedmd.com/)


#### Determining variable positions
After identifying high quality Illumina WGS data for each lab strain, we performed mapping and variant calling of each strain to the _P. falciparum_ 3D7 reference genome as described [here](https://github.com/ko-katie/pvmdr1) under "GATK Analysis", using the _P. falcpiarum_ unmasked_regions.intervals file uploaded to this repository.

The resulting combined and filtered vcf file was then used to identify positions that differentiate between the lab strains.
```
python3 parse_vcf.py /path/to/all_samples_variants_filtered2.vcf /path/to/unmasked_regions.intervals /path/to/variable_positions_output.txt /path/to/pca_matrix_output.txt
```
* Takes the combined vcf and determines positions at which not all strains share the same genotype
* all_samples_variants_filtered2.vcf is the vcf file outputted from the GATK analysis
* unmasked_regions.intervals is the path to the provided unmasked_regions.intervals file
* variable_positions_output.txt is an output file containing a list of positions at which the genotype of all of the strains are not the same
* pca_matrix_output.txt is an output file containing a matrix representation of the variable position data, used to verify that the strains can be differentiated by the identified variable positions.
```
python3 annotated_positions_by_strain.py /path/to/pca_matrix_output.txt /path/to/annotated_variable_positions_list.txt 
```
* Uses the previously outputted pca matrix to identify positions at which one strain differs from all of the others
* /path/to/pca_matrix_output.txt from running parse_vcf.py is inputted to identify positions at which the genotype of one strain differs from the genotype of all of the other strains
* /path/to/annotated_variable_positions_list.txt is the output file which contains a tab delimited file with unique positions and the strain they are unique to
```
python3 filter_35bp_window.py /path/to/annotated_variable_positions_list.txt /path/to/annotated_variable_positions_list_window_filtered.txt
```
* Takes previously generated list of positions and ensures all positions are at least 35bp apart
* /path/to/annotated_variable_positions_list.txt that was previously generated
* /path/to/annotated_variable_positions_list_window_filtered.txt output file from filtering positions to be at least 35bp apart
```
python3 get_variable_positions_and_alleles_new.py /path/to/PlasmoDB-63_Pfalciparum3D7_Genome.fasta /path/to/annotated_variable_positions_list_window_filtered.txt /path/to/all_samples_variants_filtered2.vcf /path/to/working_dir/
```
* Takes filtered list of variable positions and tests k-mers of different lengths spanning those positions to investigate necessary k-mer size for minimizing k-mers being present elsewhere in the genome
* /path/to/PlasmoDB-63_Pfalciparum3D7_Genome.fasta to determine sequences of k-mers spanning variable positions
* /path/to/annotated_variable_positions_list_window_filtered.txt is inputted to only analyze positions at least 35bp apart
* /path/to/all_samples_variants_filtered2.vcf to determine 3D7-expected and alternate alleles at each position
* /path/to/working_dir/ where two output files will be created: 1) "variable_position_alleles.txt" file storing the 3D7 and non-3D7 allele at each position and 2) "kmer_frequencies.txt" file containing the number of specific and non-specific k-mers for each strain at each k-mer length tested
```
python3 get_kmers_present_once_in_genome.py /path/to/PlasmoDB-63_Pfalciparum3D7_Genome.fasta /path/to/working_dir/variable_position_alleles.txt /path/to/25bp_kmers.txt
```
* Writes out 25bp kmers for positions where the 3D7-expected k-mer is present only once in the _P. falciparum_ genome, the non-3D7-expected k-mer is not present, and the reverse complement of either is not present.
* /path/to/PlasmoDB-63_Pfalciparum3D7_Genome.fasta to get kmers from 3D7 reference
* /path/to/working_dir/variable_position_alleles.txt to get expected specific and alternate alleles
* /path/to/25bp_kmers.txt to store correct k-mers
```
python3 filter_kmers_by_human.py /path/to/25bp_kmers.txt /path/to/25bp_kmers_human_filtered.txt /path/to/GRCh38.fasta
```
* Filter out previously generated k-mers if they or their reverse complement appear in the human genome
* /path/to/25bp_kmers.txt to filter previously generated 25bp kmers for presence in human genome
* /path/to/25bp_kmers_human_filtered.txt to output filtered 25bp kmers
* /path/to/GRCh38.fasta is used to check for k-mer presence in human genome, can be downloaded [here](https://www.gencodegenes.org/human/release_38.html)
```
python3 add_reverse_complements.py /path/to/25bp_kmers_human_filtered.txt /path/to/25bp_kmers_human_filtered_rc.txt
```
* Create final k-mer list that includes reverse complement sequence of each final k-mer
* /path/to/25bp_kmers_human_filtered.txt containing kmers not found in human
* /path/to/25bp_kmers_human_filtered_rc.txt output file with reverse complements added
