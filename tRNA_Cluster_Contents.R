# Creates a histogram that shows the number of tRNAs per cluster and whether the clusters are heterogeneous or homogeneous
# This script is tailored for Chicken

# Load packages
library(lattice)

# Import data
data<-read.table("D:/tRNAscan Analyses/tRNA Locations/clusters per species/results_cluster_analysis_chicken.txt", header = T)
head(data)

# Extract variables
scaffold<-as.factor(data$Scaffold)
locations <- c(data$begin1, data$end1)
number_of_scaffold <- nlevels(scaffold)

median_value = 129784 # Threshold for distance between clusters (calculated beforehand)
s = 1 # Number of scaffold that is being processed
all_distances <- c() # Array to save distance information
cluster_results <- {}
while (s < number_of_scaffold){
  # Extracting data per scaffold
  subset <- data[which(data$Scaffold == levels(scaffold)[s]),] # Select data with specific scaffold
  scaffold_name <- levels(scaffold)[s] # Save scaffold name
  print(scaffold_name)

  # Making new dataframe
  amino_acids <- c(as.character(subset$AA1), as.character(subset$AA2)) # Get amino acids

  begin <- c(subset$begin1, subset$begin2) # Get beginning of tRNA
  end <- c(subset$end1, subset$end2) # Get ending of tRNA

  new_data<-cbind(amino_acids, begin, end) # Combine all data into new table
  new_data<-unique(new_data) # Remove duplicates
  new_data<-new_data[order(as.numeric(new_data[,2])),] # Order by beginning of tRNAs

  if (class(new_data) == "matrix"){
    
    # Collect location data and amino acids
    locations<-sort(c(as.numeric(new_data[,2]), as.numeric(new_data[,3])))
    amino_acids_ordered <- c(new_data[,1])
    cluster<-c() # Array to save information about tRNA clusters

    x = 3 # Beginning of second AA
    y = 2 # Ending of first AA
    z = 1 # Amino Acids selection
    output <- c()
    while (z <= length(amino_acids_ordered)){
      if (z < length(amino_acids_ordered)){
        distance = locations[x] - locations[y]
        all_distances = c(all_distances, distance)
        if (distance < median_value){
          output<-c(output, amino_acids_ordered[z], "---", distance, "---") # Save distance for output-file
          cluster<-c(cluster, amino_acids_ordered[z]) # Add amino acid to vector
          } else {
          output<-c(output, amino_acids_ordered[z], "|||", distance, "|||")
          cluster<-c(cluster, amino_acids_ordered[z]) # Add amino acid to vector
          
          # Analyze cluster contents
          # Size?
          cluster_size <- length(cluster)
          # Homo- or heterogeneous content?
          if (length(unique(cluster)) == 1){
            cluster_content <- "Homogeneous"
          } else {cluster_content <- "Heterogeneous"}
          
          # Save results
          cluster_results <- rbind(cluster_results, c(scaffold_name, cluster_size, cluster_content))
          
          # Reset cluster vector
          cluster<-c()} 
        # Increment variables
        x = x + 2
        y = y + 2
        z = z + 1
      } else {
          output <- c(output, amino_acids_ordered[z])
          cluster<-c(cluster, amino_acids_ordered[z]) # Add final amino acid to vector
          cluster_size <- length(cluster)
          
          # Analyze cluster contents
          # Size
          cluster_size <- length(cluster)
          # Homo- or heterogeneous content
          if (length(unique(cluster)) == 1){
            cluster_content <- "Homogeneous"
          } else {cluster_content <- "Heterogeneous"}
          
          # Save results
          cluster_results <- rbind(cluster_results, c(scaffold_name, cluster_size, cluster_content))
          
          z = z + 1}
    }
  #sink("D:/tRNAscan Analyses/tRNA Locations/R_results_chicken.txt", append=TRUE)
  #results<-print(c(paste(scaffold_name, ":", paste(output, collapse = ""), collapse = "")), quote = FALSE)
  #sink()
  s = s + 1 # Move to next scaffold
  } else {s = s + 1}
}

# Results of cluster analysis
colnames(cluster_results) <- c("chromosome", "size", "content")
cluster_results <- cluster_results[-(1:4),] # Remove first 4 lines - scaffolds
cluster_results <- subset(cluster_results, size > 1) # Remove clusters with just one tRNA

####### Making Histogram ########
library(ggplot2)

size <- as.numeric(cluster_results[,2])
content<- as.factor(cluster_results[,3])

x<-data.frame(structure(list(size, content), .Names = c("size", "content")))

ggplot(x, aes(x=size, fill=content)) +
  stat_bin(bins = 25, binwidth = 1, col = "black") +
  scale_y_continuous(breaks=seq(0, 14, 1)) +
  scale_x_continuous(breaks=seq(1, 22, 1)) +
  labs(fill = "Cluster Content") + theme(legend.position = "bottom") +
  scale_fill_grey() +
  ggtitle("Clusters of tRNA genes in the Chicken Genome") + xlab("Number of tRNAs in Cluster") + ylab("Frequency of Cluster") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
