# Python script to collect information on genomic locations of TE-derived tRNAs
# Import files are fasta-file and text-file from BLASTn search
# Sequence similarity threshold is set by user (ranging from 1e-50 to 1e-10)

import os.path

##############################################
######## STEP 1 - GET LOCATIONS OF TE ########
#############################################

# Getting information from Fasta-file
Fasta_file = open("BLASTn_fasta.txt", "r")              # Open Fasta-file
Output_file_1 = open("output_1.txt", "a")       # Open file that will contain output: sequence_ID and scaffold name
Output_file_1.write("Sequence_ID_Fasta\tScaffold_name\n") # Header of output file 1

for line in Fasta_file:
        if ">" in line:
                contents = line.split()         # Split contents of file on whitespaces
                sequence_ID = contents[0]       # Column 0 contains the sequence ID
                sequence_ID = sequence_ID.replace(">", "") # Get rid of > sign
                if contents[5] == "chromosome":
                        scaffold_name = "".join(["chr", contents[6]])   # Column 6 contain the chromosome number
                        scaffold_name = scaffold_name.replace(",","")
                else:
                        scaffold_name = "X"

                Output_file_1.write(sequence_ID + "\t" + scaffold_name + "\n") # Write information to output-file

Output_file_1.close()
print("Collected data from Fasta-file")

# Getting information from Text-file
Text_file = open("BLASTn_text.txt", "r")
Output_file_2 = open("output_2.txt", "a")       # Open file that will contain output: sequence_ID and genomic locations
Output_file_2.write("Sequence_ID_Text\tbegin\tend\tEscore\n") # Header of output file 2

for line in Text_file:
        if "|" in line:
                contents = line.split()         # Split contents of line based on whitespace
                sequence_ID = contents[3]       # Column 3 contains sequence_ID (will be compared with Fasta-file to see if we are extracting the right information)
                location1 = contents[10]        # Column 10 contains the beginning location
                location2 = contents[11]        # Column 11 contains the ending location
                Escore = contents[12]           # Column 12 contains the E-score

                # Sort locations in begin and end
                if location1 > location2:
                        begin = location2
                        end = location1
                else:
                        begin = location1
                        end = location2

                Output_file_2.write(sequence_ID + "\t" + begin + "\t" +  end + "\t" + Escore + "\n") # Write information to output-file

Output_file_2.close()
print("Collected data from Text-file")

# Merge contents of Fasta-file and Text-file
file_1 = open("output_1.txt", "r") # Open fasta-output
file_2 = open("output_2.txt", "r") # Open text-output
TE_locations = open("output_TE_locations.txt", "a") # Open file that will contain final information

for line in file_1:
        TE_locations.write(line.rstrip() + "\t" + file_2.readline().strip() + "\n")

file_1.close()
file_2.close()
TE_locations.close()
print("Merged contents of Fasta-file and Text-file")

# Filter out sequences with low E-score
TE_locations = open("output_TE_locations.txt", "r")
Filtered_TE_locations = open("TE_locations.txt", "a")

threshold = raw_input("Which threshold would you like to use? ")
for line in TE_locations:
        contents = line.split()
        if contents[5] == "Escore":
                Filtered_TE_locations.write(line)
        elif float(contents[5]) < float(threshold):
                Filtered_TE_locations.write(line)
        else:
                next

TE_locations.close()
Filtered_TE_locations.close()
print("Sequences filtered based on E-score with a threshold of " + str(threshold))

# Remove temporary output-files
os.system("rm output_*")

# Close all files
Fasta_file.close()
Text_file.close()

######################################################
####### STEP 2 - CHECK IF tRNA OVERLAP WITH TE #######
######################################################

species = raw_input("Which species are you analyzing? ")
tRNA_file = open("private/tRNA_Overview/" + species + "_tRNA.txt", "r")
TE_tRNA = open("TE_tRNA_" + species + ".txt", "a")

for line in tRNA_file:
        # Get information on tRNAs
        contents = line.split()
        scaffold = contents[0]
        location1 = contents[2]
        location2 = contents[3]

        # Get right beginning location of tRNA
        if location1 > location2:
                begin_tRNA = location2
        else:
                begin_tRNA = location1

        # Make array of ranges for scaffold
        TE = open("TE_locations.txt", "r")
        Ranges = []
        for info in TE:
                data = info.split()
                if scaffold == data[1] :
                        TE_range = range(int(data[3]), int(data[4]))
                        Ranges.append(TE_range)
        TE.close()

        # Check if begin position of tRNA falls within one of the ranges
        for i in range(0,len(Ranges)):
                if int(begin_tRNA) in Ranges[i]:
                        TE_tRNA.write(line)
                        next

print("Filtered out TE-derived tRNAs!")

os.system("rm TE_locations.txt")

tRNA_file.close()
TE_tRNA.close()

####################################################
######### STEP 3 - SUMMARIZE RESULTS ###############
####################################################

print("\n\nSummary of TE-derived tRNAs\n")
print("Total")
os.system("cat TE_tRNA_" + species + ".txt | wc -l")

AminoAcids = ["Pseudo", "Undet", "SeC", "Sup", "Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Tyr", "Trp", "Val"]
for AA in AminoAcids:
        print(AA)
        os.system("cat TE_tRNA_" + species + ".txt | grep " + AA + "| wc -l")
