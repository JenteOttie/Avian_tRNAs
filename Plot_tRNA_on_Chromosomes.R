# R Script to plot tRNA locations on chromosomes
# https://stackoverflow.com/questions/33727432/how-to-plot-positions-along-a-chromosome-graphic

# This script is tailored for the chicken genome (other species have their own scripts)

install.packages("ggrepel")

# Load packages
library("ggplot2") # for the plot
library("ggrepel") # for spreading text labels on the plot, you can replace with `geom_text` if you want
library("scales") # for axis labels notation

##############################
######### GET DATA ###########
##############################

# Load data and create dataframe
data<-read.table("D:/tRNAscan Analyses/tRNA Locations/Clusters on Chromosomes/chicken_tRNA.txt", header = T)
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

# Putting it all together in a dataframe
sample<-data.frame(structure(list(tRNAs, chromosome, begin, end, pseudo_or_not), 
                  .Names = c("tRNA", "chromosome", "begin", "end", "pseudo")))
sample<-sample[-(286:333),] # Remove samples that do not map to a chromosome

# Get chromosome sizes (based on fasta-file galGal5.fa)
# Using this one-liner: awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' file.fa
chrom_size<-read.table("D:/tRNAscan Analyses/tRNA Locations/Clusters on Chromosomes/Chicken_Chromosomes.txt", header = T)

##############################
###### ADJUST DATA ###########
##############################

# create an ordered factor level to use for the chromosomes in all the datasets
chrom_order <- rev(c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                 "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", 
                 "chr22", "chr23", "chr24","chr25","chr26","chr27","chr28","chr30",
                 "chr31","chr32","chr33", "chrW", "chrZ"))
chrom_key <- setNames(object = as.character(c(seq(1,34))), 
                      nm = chrom_order)
chrom_order <- factor(x = chrom_order, levels = rev(chrom_order))

# convert the chromosome column in each dataset to the ordered factor
chrom_size[["chromosome"]] <- factor(x = chrom_size[["chromosome"]], 
                                      levels = chrom_order)
sample[["chromosome"]] <- factor(x = sample[["chromosome"]], 
                                     levels = chrom_order)

# create a color key for the plot
group.colors <- c(Pseudo = "red", Normal = "blue")

##############################
######## PLOT DATA ###########
##############################

ggplot(data = chrom_size) + 
  # base rectangles for the chroms, with numeric value for each chrom on the x-axis
  geom_rect(aes(xmin = as.numeric(chromosome) - 0.2, 
                xmax = as.numeric(chromosome) + 0.2, 
                ymax = size, ymin = 0), 
            colour="black", fill = "grey") + 
  # rotate the plot 90 degrees
  coord_flip() +
  # black & white color theme 
  theme(axis.text.x = element_text(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + 
  # give the appearance of a discrete axis with chrom labels
  scale_x_discrete(name = "chromosome", limits = names(chrom_key)) +
  # add lines for pseudogenes
  geom_rect(data = subset(sample, sample$pseudo == "Pseudo")
            , mapping = aes(xmin = as.numeric(chromosome) - 0.2, 
                  xmax = as.numeric(chromosome) + 0.2, 
                  ymax = end, ymin = begin, fill = pseudo), color = "red") +
  # add lines for tRNA genes
  geom_rect(data = subset(sample, sample$pseudo == "Normal")
            , mapping = aes(xmin = as.numeric(chromosome) - 0.2, 
                            xmax = as.numeric(chromosome) + 0.2, 
                            ymax = end, ymin = begin, fill = pseudo), color = "blue") +
  scale_fill_manual(name = "Genes", values = group.colors) +
  ggtitle("Chromosomal Locations of tRNA genes in the Chicken genome") +
  # supress scientific notation on the y-axis
  scale_y_continuous(labels = comma) +
  ylab("region (bp)")
