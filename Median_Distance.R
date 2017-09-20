# The median distance between consecutive tRNAs is calculated (genome-wide and chromosome-specific)

# Load data and create dataframe
data<-read.table("D:/tRNAscan Analyses/tRNA Locations/Clusters on Chromosomes/Text-files/chicken_tRNA.txt", header = T)
tRNAs<-c(as.character(data$Type))
chromosome<-c(as.character(data$Name))
begin<-c(data$Begin)
end<-c(data$End)

pseudo_or_not<-c() # Make vector to check for pseudogenes (including undetermined)
for (i in seq(1:length(tRNAs))){
  if (tRNAs[i] == "Pseudo" | tRNAs[i] == "Undet"){
    pseudo_or_not<-c(pseudo_or_not, "Pseudo")
  }
  else {pseudo_or_not<-c(pseudo_or_not, "Normal")}
}

orientation<-c() # Check for orientation of genes
for (i in seq(1:length(tRNAs))){
  if (begin[i] > end[i]){
    orientation<-c(orientation, ">")
  }
  else {orientation<-c(orientation, "<")}
}

# Putting it all together in a dataframe
all_sample<-data.frame(structure(list(tRNAs, chromosome, begin, end, pseudo_or_not, orientation), 
                             .Names = c("tRNA", "chromosome", "begin", "end", "pseudo", "orientation")))

# Filter data, get rid of pseudogenes
sample<-subset(all_sample, all_sample$pseudo == "Normal")

# Load chromosome size
chrom_size<-read.table("D:/tRNAscan Analyses/tRNA Locations/Clusters on Chromosomes/Text-files/Chicken_chromosomes.txt", header = T)

chr<-as.factor(sample$chromosome)
locations <- c(sample$begin, sample$end)
number_of_chr <- nlevels(as.factor(chromosome))

s = 1 # Number of scaffold that is being processed
all_distances <- c() # Array to save distance information
summary_results <- {} # Table to save results
while (s < number_of_chr){
  # Extracting data per scaffold
  subset <- data[which(sample$chromosome == levels(chr)[s]),] # Select data with specific scaffold
  chr_name <- levels(chr)[s] # Save scaffold name
  print(chr_name)
  
  # Extract locations and sort them
  locations <- c(subset[,3], subset[,4])
  ordered_locations<-sort(locations)
  
  if (grepl("chr", chr_name) & length(locations) > 2){ # Only use sequences that have been assigned to chromosomes and chromsomes than contain more than 1 tRNA
  
  # Calculate distances between tRNAs
  x = 3
  y = 2
  output <- c()
  chr_distances <- c() # Array to save distance information per chromosome
  while (x <= length(locations)){
    distance = ordered_locations[x] - ordered_locations[y] # Calculate distance between two tRNAs
    chr_distances = c(chr_distances, distance) # Save calculated distance in vector
    # Increment variables
    x = x + 2
    y = y + 2
  }
  s = s + 1
  print(chr_distances)
  
  num_of_tRNAs <- length(locations)/2 # Count number of tRNAs on chromosome
  med <- median(chr_distances) # Calculate median of distances
  chr_length <- chrom_size[which(chrom_size$chromosome == chr_name),2] # Get chromosome length
  summary_results <- rbind(summary_results, c(chr_name, num_of_tRNAs, med, chr_length))
  
  all_distances <- c(all_distances, chr_distances)
  }
  else {s = s + 1}
}

total_num_of_tRNA <- sum(as.numeric(summary_results[,2]))
median_total <- median(all_distances)
genome_length <- sum(chrom_size$size) # Total genome length (including chromosomes with 0 or 1 tRNAs)

summary_results <- rbind(summary_results, c("Total", total_num_of_tRNA, median_total, genome_length))
colnames(summary_results) <- c("Chromosome", "Number of tRNAs", "Median Distance", "Chromosome Length")
write.table(summary_results, "C:/Users/jente.ottenburghs/Desktop/Zebrafinch_tRNA_on_chromosomes.txt")

plot(sort(all_distances, decreasing = T))
abline(v=(length(all_distances)/2))
median(all_distances)
mean(all_distances)
