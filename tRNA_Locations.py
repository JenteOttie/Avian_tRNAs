# Python script to check the distribution of tRNAs across the genome
# The genomic locations for each tRNA are extracted from a overview file
# To check whether the tRNAs cluster in the genome, the user is asked to provide a threshold of distance between tRNAs (in this case the genome-wide median - see R-script)

import os.path

# Ask user which species they want to focus on
species = raw_input("Which species would you like to use: ")
threshold = int(raw_input("Which threshold would you like to use for the cluster analysis:  "))

# Open Overview file
input_file = open("private/tRNA_Overview/" + species + "_tRNA.txt", "r")

####################################
##### STEP 1 - Collecting data #####
####################################

# Extract locations of tRNAs per scaffold and save them in temporary files
# These temporary files will be used in the next step
number_of_tRNA = 1
scaffold_name = "X"
list_of_scaffolds = []

for line in input_file:
        if "Sequence" in line or "Name" in line or "-" in line or "Pseudo" in line or "Undet" in line: # Ignore first three lines of the tRNA file and ignore pseudogenes
                next

        else:
                line_contents = line.split('\t') # Split line according by tab
                position_1 = line_contents[2] # Beginning of tRNA
                position_2 = line_contents[3] # End of tRNA
                tRNA_type = line_contents[4] # AA coded by tRNA
                anticodon = line_contents[5] # Anticodon of tRNA

                if scaffold_name == line_contents[0]: # Compare scaffold names. If the same, collect tRNA information
                        # Collect information and save to temporary file
                        position_1 = line_contents[2] # Beginning of tRNA
                        position_2 = line_contents[3] # End of tRNA
                        tRNA_type = line_contents[4] # AA coded by tRNA
                        anticodon = line_contents[5] # Anticodon of tRNA

                        temp = open("temp_" + scaffold_name + ".txt", "a")
                        temp.write(scaffold_name + "\t" + position_1 + "\t" + position_2 + "\t" + tRNA_type + "\t" + anticodon + "\n")
                        temp.close()

                        # Counting number of tRNAs per scaffold
                        number_of_tRNA = number_of_tRNA + 1
                else:
                        #print(scaffold_name + " contains " + str(number_of_tRNA) + " tRNAs")

                        scaffold_name = line_contents[0] # Set new scaffold_name
                        list_of_scaffolds.append(scaffold_name)
                        number_of_tRNA = 1 # Reset tRNA counter

                        temp = open("temp_" + scaffold_name + ".txt", "a")
                        temp.write(scaffold_name + "\t" + position_1 + "\t" + position_2 + "\t" + tRNA_type + "\t" + anticodon + "\n") # Add first line to file
                        temp.close()

#########################################
##### STEP 2 - Check for clustering #####
#########################################

print("Data collected and stored in temporary files")
print(str(len(list_of_scaffolds)) + " scaffolds will now be analyzed")

no_tRNA = 0
begin_end_list = [] # Empty list to store begin and end locations
AA = [] # Empty list to store amino acids


# Prepare header for results file
results = open("results_cluster_analysis_" + species + ".txt", "a")
results.write("Scaffold\tDistance\tAA1\tbegin1\tend1\tAA2\tbegin2\tend2\n")
results.close()

for i in list_of_scaffolds:
        location_file = open("temp_" + i + ".txt") # Open file with tRNA locations made in Step 1
        num_lines = sum(1 for line in location_file) # Get number of lines in file
        location_file.close()

        # Remove files with just one tRNA
        if (num_lines == 1): # If the file contains only one tRNA discard it.
                print("scaffold " + i + " contains only one tRNA and will be removed")
                os.remove("temp_" + i + ".txt") # Remove file
                no_tRNA = no_tRNA + 1 # Keep track of files with only one tRNA
        else:
                print("Analyzing scaffold " + i + ", which contains " + str(num_lines) + " tRNAs" )
                location_file = open("temp_" + i + ".txt") # Re-open file
                for data in location_file:
                        data_contents = data.split('\t') # Split line according by tab
                        begin = int(data_contents[1]) # Beginning of tRNA
                        end = int(data_contents[2]) # End of tRNA
                        amino = data_contents[3] # Amino Acid
                        anticodon = data_contents[4].rstrip("\n") # Anticodon with function rstrip() to remove newline (\n) at the end

                        begin_end_list.append(begin) # Add begin and end locations to list
                        begin_end_list.append(end)
                        AA.append(amino) # Save amino acids in list AA. Do it twice to make extracting easier in later stage
                        AA.append(amino)

        # Analyzing files with more than one tRNA
        if len(begin_end_list) == 0:
                next
        else:
                [x for y,x in sorted(zip(AA,begin_end_list))]   # Sort lists for amino acids and locations at the same time
                                                                # They have to be sorted together, otherwise the order of the AAs gets lost!

                # Calculate distance between end of one tRNA and beginning of the next one
                x = 1
                y = 2
                while y < len(begin_end_list):
                        distance = begin_end_list[y] - begin_end_list[x]
                        if (distance < threshold and distance > -threshold): # Save tRNA genes that are less than 1000 bp apart
                                results = open("results_cluster_analysis_" + species + ".txt", "a")
                                results.write(i + "\t" + str(distance) + "\t" + AA[x] + "\t" + str(begin_end_list[x-1]) + "\t" + str(begin_end_list[x]) + "\t" + AA[y] + "\t" + str(begin_end_list[y]) + "\t" + str(begin_end_list[y+1]) + "\n")
                                results.close()
                        else:
                                next
                        y  = y + 2
                        x  = x + 2

        begin_end_list = [] # Reset lists for next scaffold
        AA = []

# Delete temporary files
os.system("rm temp_*")
