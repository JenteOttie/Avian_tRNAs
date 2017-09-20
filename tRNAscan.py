# The bird genomes are spread over three different folders on the UPPMAX server
# 48_genomes: containing the 48 genomes used in the Science paper by Jarvis et al. (2014)
# more_bird_genomes: containing genomes published after the Science paper but before June 2016 (added by Alexander Suh)
# extra_genomes: containing genomes published after June 2016 (added by Jente Ottenburghs)

# This Python script runs tRNAscan-SE on all the available bird genomes

################################################################
############## run tRNAscan on 48_genomes folder ###############
################################################################

import os.path

# List of species to be analyzed
species = ["white-tailed_eagle", "bald_eagle", "vulture", "duck", "swift", "hummingbird", "hornbill", "seriema", "killdeer", "mousebird", "pigeon", "bee-eater", "cuckoo", "sunbittern", "peregrine", "loon", "crane", "cuckoo-roller", "mesite", "turaco", "hoatzin", "bustard", "rifleman", "american_crow", "zebrafinch", "golden-collared_manakin", "groundfinch", "egre", "pelican", "cormorant", "ibis", "tropicbird", "flamingo", "woodpecker", "grebe", "fulmar", "kea", "sandgrouse", "emperor_penguin", "adelie_penguin", "owl", "tinamou", "trogon"]

# Names of files in the folder (corresponding to the species list above)
code = ["hala", "hall", "cath", "anas", "chae", "caly", "buce", "cari", "char", "coli", "colu", "mero", "cucu", "eury", "falc", "gavi", "bale", "lept", "mesi", "taur", "opis", "chla", "acan", "corv", "taen", "mana", "geos", "egre", "pele", "phal", "nipp", "phae", "phoe", "pico", "podi", "fulm", "nest", "pter", "apte", "pygo", "tyto", "tina", "apal"]

for i in range(0, len(species)):
	print( "Analyzing " + species[i] + " with code " + code[i] + "\n")

	#Unzip Fastafile
	os.system("gunzip /proj/snic2017-7-108/private/48_genomes/" + code[i] + ".fa.gz")
	#Run tRNAscan
	# -o saves an overview file with genomic locations of all tRNAs
	# -m saves a summary file on the tRNA content
	# -f saves a file containing sequence information of each tRNA
	os.system("tRNAscan-SE -o /home/jente/private/tRNA_Overview/" + species[i] + "_tRNA.txt -m /home/jente/private/tRNA_Summaries" + species[i] + ".txt -f /home/jente/private/tRNA_Sequences/" + species[i] + "_sequences.txt  /proj/snic2017-7-108/private/48_genomes/" + code[i] + ".fa")
	#Re-zip the Fasta-file
	os.system("gzip /proj/snic2017-7-108/private/48_genomes/" + code[i] + ".fa")
  
################################################################################
############## run tRNAscan on more_bird_genomes_June2016 folder ###############
################################################################################

import os.path

# List of species to be analyzed
species = ["eagle", "goose", "kiwi", "ruff", "saker_falcon", "bobwhite", "crow", "white-throated_sparrow", "canary", "flycatcher", "great_tit", "blue_tit", "ground_tit", "plumbeitarsus", "trochiloides", "viridanus", "blue-crowned_manakin", "starling", "silveye", "blue_amazon", "puerto_rico", "macaw", "ostrich"]

# Names of files in the folder (corresponding to the species list above)
code = ["aquChr1", "ansCyg1", "aptMan1", "calPug_ruff_EG1", "falChe1", "colVir1", "corCor1", "zonAlb1", "serCan1", "ficAlb2", "parMaj1", "cyaCae2", "pseHum1", "phyPlu1", "phyTro1", "phyVir1", "lepCor1", "stuVul1", "zosLat1", "amaAes1", "amaVit1", "araMac1", "strCam2"]

for i in range(0, len(species)):
	print( "Analyzing " + species[i] + " with code " + code[i] + "\n")

	#Unzip Fastafile
	os.system("gunzip /proj/snic2017-7-108/private/more_bird_genomes_June2016/" + code[i] + ".fa.gz")
	#Run tRNAscan
        # -o saves an overview file with genomic locations of all tRNAs
        # -m saves a summary file on the tRNA content
        # -f saves a file containing sequence information of each tRNA
        os.system("tRNAscan-SE -o /home/jente/private/tRNA_Overview/" + species[i] + "_tRNA.txt -m /home/jente/private/tRNA_Summaries" + species[i] + ".txt -f /home/jente/private/tRNA_Sequences/" + species[i] + "_sequences.txt  /proj/snic2017-7-108/private/more_bird_genomes_June2016/" + code[i] + ".fa")
	#Re-zip the Fasta-file
	os.system("gzip /proj/snic2017-7-108/private/more_bird_genomes_June2016/" + code[i] + ".fa")
  
###################################################################
############## run tRNAscan on extra_genomes folder ###############
###################################################################

import os.path

# List of all the species that will be analyzed
species = ["spot-billed_duck", "nightjar", "stork", "band-tailed_pigeon", "guineafowl", "scaled_quail", "turkey", "grus", "munia", "warbler", "sparrow", "seedeater", "auritus", "brasilianus", "harrisi", "urile", "budgerigar"]

# Names of files in the folder (corresponding to the species list above)
code = ["anas_zonorhyncha", "antrostomus", "ciconia", "patagoenias", "numida", "callipepla", "meleagris", "grus", "lonchura", "setophaga", "passer", "sporophila", "nannopterum_auritus", "nannopterum_brasilianus", "nannopterum_harrisi", "urile", "melopsittacus"]

for i in range(0, len(species)):
	print( "Analyzing " + species[i] + " with code " + code[i] + "\n")

	#Unzip Fastafile
	os.system("gunzip /proj/snic2017-7-108/private/extra_bird_genomes/" + code[i] + ".fa.gz")
	#Run tRNAscan
        # -o saves an overview file with genomic locations of all tRNAs
        # -m saves a summary file on the tRNA content
        # -f saves a file containing sequence information of each tRNA
        os.system("tRNAscan-SE -o /home/jente/private/tRNA_Overview/" + species[i] + "_tRNA.txt -m /home/jente/private/tRNA_Summaries" + species[i] + ".txt -f /home/jente/private/tRNA_Sequences/" + species[i] + "_sequences.txt  /proj/snic2017-7-108/private/extra_bird_genomes/" + code[i] + ".fa")
	#Re-zip the Fasta-file
	os.system("gzip /proj/snic2017-7-108/private/extra_bird_genomes/" + code[i] + ".fa")
