# Heatmap of tRNA content on amino acid and anticodon levels
#https://www.r-bloggers.com/heatmaply-interactive-heat-maps/

################################################
### Installing and loading required packages ###
################################################
install.packages("heatmaply")
library(ggplot2)
library(plotly)
library(heatmaply)
library(openxlsx)

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}


##################### AMINO ACIDS #########################

###########################################################
### Reading in data and transform it into matrix format ###
###########################################################

library(xlsx)

data <- read.xlsx("D:/tRNAscan Analyses/Overview.xlsx", sheetName = 'Amino Acids')
rnames <- data[,3]                              # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,7:27]) # transform column 2-n into a matrix
rownames(mat_data) <- rnames                    # assign row names

#############################################################
### Customizing and plotting the heat map for Amino Acids ###
#############################################################

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("lightyellow", "lightblue", "red"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0,33,length=100),   # for yellow
               seq(34,66,length=100),  # for blue
               seq(67,100,length=100)) # for red

# creates a 5 x 5 inch image
png("heatmaps_in_r.png",    # create PNG for the heat map
    width = 5*600,        # 5 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 5)        # smaller font size

# Create heatmap that keeps the original order of the data
heatmap.2(mat_data,
          #cellnote = mat_data, # same data set for cell labels
          main = "Heatmap",     # heat map title
          notecol="black",      # change font color of cell labels to black
          cexRow = 0.5,         # Font size
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA",            # turn off column clustering
          Rowv = FALSE)         # keep original order of data

dev.off() 

# Create heatmap that clusters species according to AA content (great to see outliers!)
heatmap.2(mat_data,
          #cellnote = mat_data, # same data set for cell labels
          main = "Heatmap",     # heat map title
          notecol="black",      # change font color of cell labels to black
          cexRow = 0.5,         # Font size
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA")            # turn off column clustering         
#Colv = as.dendrogram(cluster)) # apply default clustering method
dev.off() 


##################### ANTICODONS ##########################

###########################################################
### Reading in data and transform it into matrix format ###
###########################################################

library(xlsx)

data <- read.xlsx("D:/tRNAscan Analyses/Overview.xlsx", sheetName = 'No AA Header', header = T)
rnames <- data[,3]                              # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,7:68]) # transform column 2-n into a matrix
rownames(mat_data) <- rnames                    # assign row names

#############################################################
### Customizing and plotting the heat map for Amino Acids ###
#############################################################

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("lightyellow", "lightblue", "red"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0,33,length=100),   # for yellow
               seq(34,66,length=100),  # for blue
               seq(67,100,length=100)) # for red

# creates a 5 x 5 inch image
png("heatmaps_in_r.png",    # create PNG for the heat map
    width = 5*600,        # 5 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 5)        # smaller font size

# Create heatmap that keeps the original order of the data
heatmap.2(mat_data,
          #cellnote = mat_data, # same data set for cell labels
          main = "Heatmap",     # heat map title
          notecol="black",      # change font color of cell labels to black
          cexRow = 0.5,         # Font size
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA",            # turn off column clustering
          Rowv = FALSE)         # keep original order of data

dev.off()

# Create heatmap that clusters species according to AA content (great to see outliers!)
heatmap.2(mat_data,
          #cellnote = mat_data, # same data set for cell labels
          main = "Heatmap",     # heat map title
          notecol="black",      # change font color of cell labels to black
          cexRow = 0.5,         # Font size
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA")            # turn off column clustering         

dev.off() 
