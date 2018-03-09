###### This script is used to build an input for the Clustering analysis #########
### Author: Meghna Verma
### Date: October 14, 2017, NIMMLab; Modified: Nov 7, 2017
## Last modified: January 9, 2018; Feb 12, 2018, With inclusion for data regarding patterns
## Last modified: March 9, 2017; New data files for analysis

### Input: All csv files from SFX, average values
### Output: Input for clustering, cluster numbers; Regression based analysis

# Working directory -> Set ot File panel location
library("qpcR")
#### Reading all csv files ####
temp = list.files(pattern="*.csv")

### Read 1st csv file ###
Th17 = read.csv(temp[1], header = FALSE)
Th17 = Th17[,2];
### Append the rest to it ###
for (i in 2:length(temp))
{
  file_one <- read.csv(temp[i],header = FALSE)
  file_one <- file_one[,2]
  ### This pads NA to the rows that are shorter than rest ###
  Th17 <- qpcR:::cbind.na(Th17, file_one)
}

df_Th17 <- data.frame(Th17)

write.csv(df_Th17, "Th17_LP_input.csv")


### Option 1 ####
library(e1071)

##Specific to HPLP##
#input_tobeclustered2<-read.csv("/Users/meghnaverma/Documents/MATLAB/NIMML/ENISI_output_SA_avg/LP_HP/cluser_files/LP_HP_input.csv", sep=",", header = TRUE)
#input_tobeclustered2 <- read.csv("/Users/meghnaverma/Documents/MATLAB/NIMML/Round2_SA_Xi/ENISI_Avg/Tr_LP/Cluster_files/Tr_LP_input.csv", sep=",", header = TRUE)
input_tobeclustered2 <- read.csv("/Users/meghnaverma/Documents/MATLAB/NIMML/SA-completeanalysis/ENISI-Avg/LP-Th17/CLusterfiles/Th17_LP_input.csv", sep=",", header = TRUE)

input_tobeclustered2 <- input_tobeclustered2[,-1] # remove 1:35 when you want all the time points

input_tobeclustered <- t(input_tobeclustered2)
transposed <- input_tobeclustered
#transposed <- t(input_tobeclustered)
transposed[is.na(transposed)] <- 0
#transposed <- transposed[-1,]

###CMEANS Clustering here #####
clustered_data<-cmeans(transposed, 3, iter.max = 100, verbose = FALSE, dist = "euclidean",
                       method = "cmeans", m = 2, rate.par = NULL, weights = 1)
# 
write.csv(clustered_data$cluster, "Th17_LP_clusters_numbers.csv")
write.csv(clustered_data$membership, "Th17_LPclusters.csv")
# x11()
# library(corrplot)
# corrplot(t(clustered_data$membership[1:30,]), is.corr = FALSE)

#### HCLUST HERE ####
library(rafalib)
hc <- hclust(dist(input_tobeclustered)^2, "cen")
# define some clusters

mycl <- cutree(hc, k = 3)

# get a color palette equal to the number of clusters
clusterCols <- rainbow(length(unique(mycl)))

# create vector of colors for side bar
myClusterSideBar <- clusterCols[mycl]

# choose a color palette for the heat map
myheatcol <- rev(redgreen(35))

###
#if(!require(dendextend)) install.packages("dendextend")
library(dendextend)

par(mfrow = c(1,2), mar = c(5,2,1,0))

dend <- as.dendrogram(hc)

dend <- dend %>%
  color_branches(k = 3) %>%
  set("branches_lwd", c(2,1,2)) %>%
  set("branches_lty", c(1,2,1))

#plot(dend)

dend <- color_labels(dend, k = 3)
# The same as:
# labels_colors(dend)  <- get_leaves_branches_col(dend)
#plot(dend)

col_labels <- get_leaves_branches_col(dend)

col_labels <- col_labels[order(order.dendrogram(dend))]

library(gplots)
#color.palette  <- colorRampPalette(c("#FC8D59", "#FFFFBF", "#91CF60"))(length(palette.breaks) - 1)
#my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
#my_palette <- colorRampPalette(c('red','yellow','green'))(256)
##
# draw the heat map

library("RColorBrewer")
#display all colour schemes
display.brewer.all()
colfunc <- colorRampPalette(c("#FC8D59", "#FFFFBF", "#91CF60"))
#brewer.pal(11,"RdBu")
x11()
heatmap.2(input_tobeclustered, main="Hierarchical Cluster", Rowv=dend, Colv=NA, dendrogram="row", scale="row", col=colfunc(45), density.info="none", trace="none", RowSideColors= col_labels, # to add nice colored strips        
          colRow = col_labels)
#### Making clusters for HPLP #####
## Gets all the inputs from clustered outputs and plots clusters on same plots##
#Author: Meghna Verma
#Date: July 26, 2017 (Last modified)
#NIMML ENISI MSM SA analysis - PLotting the clusters of designs for each output 

## Input file: 1. Cluster csv (number) 'likelihood' file obtained from Clustering code
##             2. Parameter design file
##             3. All .csv files with avg values (in the current directory)
##             4. Corresponding cells file 

## Output file: 1. Csv files that are clustered in different sets
##              2. Graphs that show the different clusters     

# Import the file
cluster1 = vector(); # creates an empty variable to be filled up
cluster2 = vector();
cluster3 = vector();
# 
# #install.packages("readr")
library(readr)
#dataset <- read_csv("/Users/meghnaverma/Documents/MATLAB/NIMML/Round2_SA_Xi/ENISI_Avg/iTreg_LP/Cluster_files/iTreg_LP_clusters_numbers.csv", 
                    # col_types = cols(X1 = col_character()))

dataset <- read_csv("/Users/meghnaverma/Documents/MATLAB/NIMML/SA-completeanalysis/ENISI-Avg/LP-Tr/Clusterfiles/Tr_LP_clusters_numbers.csv", col_types = cols(X1 = col_character()))

# 
#Parameters_old <- read_csv("/Users/meghnaverma/Documents/MATLAB/NIMML/ENISI_output_SA_avg/LP_Macrophage_State2/cluster_numbers/Parameters.csv", col_names = FALSE)
Parameters <- read_csv("/Users/meghnaverma/Documents/MATLAB/NIMML/SA-completeanalysis/Parameters-Combined-Regression.csv",col_names = FALSE)

# 
sorted <- dataset[with(dataset, order(x)),]
# 
index <- sort.int(dataset$x, index.return = TRUE)
# 
for (i in 1:267)
{
   if (sorted$x[i] == 1)
   {
     cluster1 <- rbind (cluster1,index$ix[i])
   }
   else if (sorted$x[i] == 2)
   {
     cluster2 <- rbind (cluster2,index$ix[i])
   }
   else
   {
     cluster3 <- rbind (cluster3,index$ix[i])
   }
 }

# # Make vectors of names of csv files that correspond to the clustered sets.
 cluster1_files <- vector()
 cluster2_files <- vector()
 cluster3_files <- vector()
# 
para1_files <- vector()
para2_files <- vector()
para3_files <- vector()
 
temp = list.files(pattern="*.csv")
#Correponding_cells <- read_csv("/Users/meghnaverma/Documents/MATLAB/NIMML/ENISI_output_SA_avg/LP_Macrophage_State2/cluster_numbers/Correponding_cells.csv")
 
#Correponding_cells <- read_csv("/Users/meghnaverma/Documents/MATLAB/NIMML/SA-completeanalysis/Correponding_cells-NEW.csv")
Correponding_cells <- read_csv("/Users/meghnaverma/Documents/MATLAB/NIMML/SA-completeanalysis/Correponding_cells-March8.csv")


SAMPLE1 <- Correponding_cells
sorted_matlab <- sort.int(SAMPLE1$MATLAB, index.return = TRUE)

for (j in 1:length(cluster1))
 {
   #print (d[j])
   cluster1_files <- rbind (cluster1_files, temp[cluster1[j]])
   para1_files <- rbind(para1_files, Parameters[sorted_matlab$ix[cluster1[j]],])
 }
for (k in 1:length(cluster2))
 {
  #print (d[j])
  cluster2_files <- rbind (cluster2_files, temp[cluster2[k]])
   para2_files <- rbind(para2_files, Parameters[sorted_matlab$ix[cluster2[k]],])
 }
 
for (l in 1:length(cluster3))
 {
#   #print (d[j])
   cluster3_files <- rbind (cluster3_files, temp[cluster3[l]])
   para3_files <- rbind(para3_files, Parameters[sorted_matlab$ix[cluster3[l]],]) 
 }
# 
# #######
# 
para1_patterns <- para1_files
para2_patterns <- para2_files
para3_patterns <- para3_files


# #files <- list.files(pattern = ".csv") ## creates a vector with all file names in your folder
# 
# ##Cluster 1 mean
polmean1 <- rep(0,length(cluster1_files))
 for(i in 1:length(cluster1_files)){
   data_my1 <- read.csv(cluster1_files[i],header=FALSE)
   data_my1 <- data_my1[,2]
   polmean1[i] <- mean(data_my1)
 }
# ##Cluster 2 mean
 polmean2 <- rep(0,length(cluster2_files))
 for(i in 1:length(cluster2_files)){
   data_my2 <- read.csv(cluster2_files[i],header=FALSE)
   data_my2 <- data_my2[,2]
   polmean2[i] <- mean(data_my2)
 }
# 
# ##Cluster 3 mean
 polmean3 <- rep(0,length(cluster3_files))
 for(i in 1:length(cluster3_files))
   {
   data_my3 <- read.csv(cluster3_files[i],header=FALSE)
   data_my3 <- data_my3[,2]
   polmean3[i] <- mean(data_my3)
 }
# 
 para1_files$mean <- polmean1
 para2_files$mean <- polmean2
 para3_files$mean <- polmean3
# 
 x <- data.frame(para1_files)
 y <- data.frame(para2_files)
 z <- data.frame(para3_files)
# 
 indexes_info <- rbind(cluster1, cluster2, cluster3)
# 
 final <- rbind(x,y,z)
 
 ### BEGIN for Pattern in last column ###
 ## Used for putting pattern as last column ##
 list1 <- 1;
 list2 <- 2;
 list3 <- 3;
 
 pattern1 <- rep("1", nrow(para1_patterns))
 para1_patterns$pattern <- pattern1
 pattern2 <- rep("2", nrow(para2_patterns))
 para2_patterns$pattern <- pattern2
 pattern3 <- rep("3", nrow(para3_patterns))
 para3_patterns$pattern <- pattern3
 
 x1_pattern <- data.frame(para1_patterns)
 y1_pattern <- data.frame(para2_patterns)
 z1_pattern <- data.frame(para3_patterns)
 
 final1_pattern <- rbind(x1_pattern,y1_pattern,z1_pattern)
 final1_pattern <- final1_pattern[,-1]
 
 #write.csv(MyData, file = "MyData.csv",row.names=FALSE)
 write.csv(final1_pattern, file = 'LP_Tr_final.csv', row.names = FALSE)
### END of pattern ####
  #####
 
 #final$index <- indexes_info
# 
# #write.csv(MyData, file = "MyData.csv",row.names=FALSE)
#write.csv(final, file = '/Users/meghnaverma/Documents/MATLAB/NIMML/Round2_SA_Xi/ENISI_Avg/iTreg_LP/iTreg_LP_means.csv', row.names = FALSE)
write.csv(final, file = '/Users/meghnaverma/Documents/MATLAB/NIMML/SA-completeanalysis/ENISI-Avg/LP-Th17/LP_Th17_means.csv', row.names = FALSE)

#####CODE FOR PLOT CLUSTERS ENDS #######

##### CODE FOR PRCC BEGINS ######

##Now normalize the data and run PCC
#data_LP_Mreg_main <- read.csv('/Users/meghnaverma/Documents/MATLAB/NIMML/ENISI_output_SA_avg/Lumen_HP/Lumen_HP_means.csv', header = TRUE)

data_LP_HP_main <- read.csv('/Users/meghnaverma/Documents/MATLAB/NIMML/SA-completeanalysis/ENISI-Avg/LP-HP/LP_HP_means.csv', header = TRUE)

data_LP_Mres_main <- read.csv('/Users/meghnaverma/Documents/MATLAB/NIMML/SA-completeanalysis/ENISI-Avg/LP-Mres/LP_Mres_means.csv', header = TRUE)
data_LP_Mreg_main <- read.csv('/Users/meghnaverma/Documents/MATLAB/NIMML/SA-completeanalysis/ENISI-Avg/LP-Mreg/LP_Mreg_means.csv', header = TRUE)
data_LP_Minf_main <- read.csv('/Users/meghnaverma/Documents/MATLAB/NIMML/SA-completeanalysis/ENISI-Avg/LP-Minf/LP_Minf_means.csv', header = TRUE)

data_LP_Th1_main <- read.csv('/Users/meghnaverma/Documents/MATLAB/NIMML/SA-completeanalysis/ENISI-Avg/LP-Th1/LP_Th1_means.csv', header = TRUE)
data_LP_Th17_main <- read.csv('/Users/meghnaverma/Documents/MATLAB/NIMML/SA-completeanalysis/ENISI-Avg/LP-Th17/LP_Th17_means.csv', header = TRUE)
data_LP_iTreg_main <- read.csv('/Users/meghnaverma/Documents/MATLAB/NIMML/SA-completeanalysis/ENISI-Avg/LP-iTreg/LP_iTreg_means.csv', header = TRUE)
data_LP_Tr_main <- read.csv('/Users/meghnaverma/Documents/MATLAB/NIMML/SA-completeanalysis/ENISI-Avg/LP-Tr/LP_Tr_means.csv', header = TRUE)

data_GLN_Th1_main <- read.csv('/Users/meghnaverma/Documents/MATLAB/NIMML/SA-completeanalysis/ENISI-Avg/GLN-Th1/GLN_Th1_means.csv', header = TRUE)
data_GLN_Th17_main <- read.csv('/Users/meghnaverma/Documents/MATLAB/NIMML/SA-completeanalysis/ENISI-Avg/GLN-Th17/GLN_Th17_means.csv', header = TRUE)
data_GLN_iTreg_main <- read.csv('/Users/meghnaverma/Documents/MATLAB/NIMML/SA-completeanalysis/ENISI-Avg/GLN-iTreg/GLN_iTreg_means.csv', header = TRUE)

data_GLN_eDC_main <- read.csv('/Users/meghnaverma/Documents/MATLAB/NIMML/SA-completeanalysis/ENISI-Avg/GLN-eDC/GLN_eDC_means.csv', header = TRUE)
#data_GLN_tDC_main <- read.csv('/Users/meghnaverma/Documents/MATLAB/NIMML/Round2_SA_Xi/ENISI_Avg/tDC_GLN/GLN_tDC_means.csv', header = TRUE)

data_GLN_tDC_main <- read.csv('/Users/meghnaverma/Documents/MATLAB/NIMML/SA-completeanalysis/ENISI-Avg/GLN-tDc/GLN_tDC_means.csv', header = TRUE)


## Remove the mean and 1st column with indices of parameters 
sorted_data_LPHP <- data_LP_HP_main[order(data_LP_HP_main$X1),] 
sorted_data_LP_Mres <- data_LP_Mres_main[order(data_LP_Mres_main$X1),] 
sorted_data_LPMreg <- data_LP_Mreg_main[order(data_LP_Mreg_main$X1),] 
sorted_data_LPMinf <- data_LP_Minf_main[order(data_LP_Minf_main$X1),] 

sorted_data_LPTh1 <- data_LP_Th1_main[order(data_LP_Th1_main$X1),] 
sorted_data_LPTh17 <- data_LP_Th17_main[order(data_LP_Th17_main$X1),] 
sorted_data_iTreg <- data_LP_iTreg_main[order(data_LP_iTreg_main$X1),] 
sorted_data_Tr <- data_LP_Tr_main[order(data_LP_Tr_main$X1),] 

sorted_data_Th1g <- data_GLN_Th1_main[order(data_GLN_Th1_main$X1),] 
sorted_data_Th17g <- data_GLN_Th17_main[order(data_GLN_Th17_main$X1),] 
sorted_data_iTregg <-data_GLN_iTreg_main[order(data_GLN_iTreg_main$X1),] 

sorted_data_eDC <-data_GLN_eDC_main[order(data_GLN_eDC_main$X1),] 
sorted_data_tDC <- data_GLN_tDC_main[order(data_GLN_tDC_main$X1),] 

#### 
data_LP_HP <- data.frame(sorted_data_LPHP[,-40])
#data_LP_HP <- data.frame(data_LP_HP[,-20])
data_LP_HP <- data.frame(data_LP_HP[,-1])

data_LP_Mres <- data.frame(sorted_data_LP_Mres[,-40]) 
#data_LP_Mres <- data.frame(data_LP_Mres[,-20])
data_LP_Mres <- data.frame(data_LP_Mres[,-1])

data_LP_Mreg <- data.frame(sorted_data_LPMreg[,-40]) 
#data_LP_Mreg <- data.frame(data_LP_Mreg[,-20])
data_LP_Mreg <- data.frame(data_LP_Mreg[,-1])

data_LP_Minf <- data.frame(sorted_data_LPMinf[,-40]) 
#data_LP_Minf <- data.frame(data_LP_Minf[,-20])
data_LP_Minf <- data.frame(data_LP_Minf[,-1])

data_LP_Th1 <- data.frame(sorted_data_LPTh1[,-40]) 
#data_LP_Th1 <- data.frame(data_LP_Th1[,-20])
data_LP_Th1 <- data.frame(data_LP_Th1[,-1])

data_LP_Th17 <- data.frame(sorted_data_LPTh17[,-40]) 
#data_LP_Th17 <- data.frame(data_LP_Th17[,-20])
data_LP_Th17 <- data.frame(data_LP_Th17[,-1])

data_LP_iTreg <- data.frame(sorted_data_iTreg[,-40]) 
#data_LP_iTreg <- data.frame(data_LP_iTreg[,-20])
data_LP_iTreg <- data.frame(data_LP_iTreg[,-1])

data_LP_Tr <- data.frame(sorted_data_Tr[,-40]) 
#data_LP_Tr <- data.frame(data_LP_Tr[,-20])
data_LP_Tr <- data.frame(data_LP_Tr[,-1])

data_GLN_Th1 <- data.frame(sorted_data_Th1g[,-40]) 
#data_GLN_Th1 <- data.frame(data_GLN_Th1[,-20])
data_GLN_Th1 <- data.frame(data_GLN_Th1[,-1])

data_GLN_Th17 <- data.frame(sorted_data_Th17g[,-40]) 
#data_GLN_Th17 <- data.frame(data_GLN_Th17[,-20])
data_GLN_Th17 <- data.frame(data_GLN_Th17[,-1])

data_GLN_iTreg <- data.frame(sorted_data_iTregg[,-40]) 
#data_GLN_iTreg <- data.frame(data_GLN_iTreg[,-20])
data_GLN_iTreg <- data.frame(data_GLN_iTreg[,-1])

data_GLN_eDC <- data.frame(sorted_data_eDC[,-40]) 
#data_GLN_eDC <- data.frame(data_GLN_eDC[,-20])
data_GLN_eDC <- data.frame(data_GLN_eDC[,-1])

data_GLN_tDC <- data.frame(sorted_data_tDC[,-40]) 
#data_GLN_tDC <- data.frame(data_GLN_tDC[,-20])
data_GLN_tDC <- data.frame(data_GLN_tDC[,-1])


## normalize data in columns 
scaled.dat_HP_LP <- scale(data_LP_HP)
colMeans(scaled.dat_HP_LP)
apply(scaled.dat_HP_LP, 2, sd)

scaled.dat_Mres <- scale(data_LP_Mres)
colMeans(scaled.dat_Mres)
apply(scaled.dat_Mres, 2, sd)

scaled.dat_Mreg <- scale(data_LP_Mreg)
colMeans(scaled.dat_Mreg)
apply(scaled.dat_Mreg, 2, sd)

scaled.dat_Minf <- scale(data_LP_Minf)
colMeans(scaled.dat_Minf)
apply(scaled.dat_Minf, 2, sd)

scaled.dat_Th1lp <- scale(data_LP_Th1)
colMeans(scaled.dat_Th1lp)
apply(scaled.dat_Th1lp, 2, sd)

scaled.dat_Th17lp <- scale(data_LP_Th17)
colMeans(scaled.dat_Th17lp)
apply(scaled.dat_Th17lp, 2, sd)

scaled.dat_iTreglp <- scale(data_LP_iTreg)
colMeans(scaled.dat_iTreglp)
apply(scaled.dat_iTreglp, 2, sd)

scaled.dat_Trlp <- scale(data_LP_Tr)
colMeans(scaled.dat_Trlp)
apply(scaled.dat_Trlp, 2, sd)

scaled.dat_Th1gln <- scale(data_GLN_Th1)
colMeans(scaled.dat_Th1gln)
apply(scaled.dat_Th1gln, 2, sd)

scaled.dat_Th17gln <- scale(data_GLN_Th17)
colMeans(scaled.dat_Th17gln)
apply(scaled.dat_Th17gln, 2, sd)

scaled.dat_iTreggln <- scale(data_GLN_iTreg)
colMeans(scaled.dat_iTreggln)
apply(scaled.dat_iTreggln, 2, sd)


scaled.dat_eDC <- scale(data_GLN_eDC)
colMeans(scaled.dat_eDC)
apply(scaled.dat_eDC, 2, sd)

scaled.dat_tDC <- scale(data_GLN_tDC)
colMeans(scaled.dat_tDC)
apply(scaled.dat_tDC, 2, sd)

#####
total_data_HPLP<- data.frame(cbind(scaled.dat_HP_LP, sorted_data_LPHP[,40]))
total_data_Mres <- data.frame(cbind(scaled.dat_Mres, sorted_data_LP_Mres[,40]))
total_data_Mreg <- data.frame(cbind(scaled.dat_Mreg, sorted_data_LPMreg[,40]))
total_data_Minf <- data.frame(cbind(scaled.dat_Minf, sorted_data_LPMinf[,40]))
total_data_Th1 <- data.frame(cbind(scaled.dat_Th1lp, sorted_data_LPTh1[,40]))
total_data_Th17 <- data.frame(cbind(scaled.dat_Th17lp, sorted_data_LPTh17[,40]))
total_data_iTreg <- data.frame(cbind(scaled.dat_iTreglp, sorted_data_iTreg[,40]))
total_data_Tr<- data.frame(cbind(scaled.dat_Trlp, sorted_data_Tr[,40]))
total_data_Th1g <- data.frame(cbind(scaled.dat_Th1gln, sorted_data_Th1g[,40]))
total_data_Th17g <- data.frame(cbind(scaled.dat_Th17gln, sorted_data_Th17g[,40]))
total_data_iTregg <- data.frame(cbind(scaled.dat_iTreggln, sorted_data_iTregg[,40]))
total_data_eDC <- data.frame(cbind(scaled.dat_eDC, sorted_data_eDC[,40]))
total_data_tDC <- data.frame(cbind(scaled.dat_tDC, sorted_data_tDC[,40]))

library("sensitivity")
coeff_pcc_HPLP <- pcc(total_data_HPLP[,1:38], total_data_HPLP[,39], nboot = 500, rank = TRUE, conf = 0.95)
coeff_pcc_Mres <- pcc(total_data_Mres[,1:38], total_data_Mres[,39], nboot = 500, rank = TRUE, conf = 0.95)
coeff_pcc_Mreg <- pcc(total_data_Mreg[,1:38], total_data_Mreg[,39], nboot = 500, rank = TRUE, conf = 0.95)
coeff_pcc_Minf <- pcc(total_data_Minf[,1:38], total_data_Minf[,39], nboot = 500, rank = TRUE, conf = 0.95)

coeff_pcc_Th1 <- pcc(total_data_Th1[,1:38], total_data_Th1[,39], nboot = 500, rank = TRUE, conf = 0.95)
coeff_pcc_Th17 <- pcc(total_data_Th17[,1:38], total_data_Th17[,39], nboot = 500, rank = TRUE, conf = 0.95)
coeff_pcc_iTreg <- pcc(total_data_iTreg[,1:38], total_data_iTreg[,39], nboot = 500, rank = TRUE, conf = 0.95)
coeff_pcc_Tr <- pcc(total_data_Tr[,1:38], total_data_Tr[,39], nboot = 300, rank = TRUE, conf = 0.95)

coeff_pcc_eDC <- pcc(total_data_eDC[,1:38], total_data_eDC[,39], nboot = 500, rank = TRUE, conf = 0.95)
coeff_pcc_tDC <- pcc(total_data_tDC[,1:38], total_data_tDC[,39], nboot = 500, rank = TRUE, conf = 0.95)


coeff_pcc_Th1g <- pcc(total_data_Th1g[,1:38], total_data_Th1g[,39], nboot = 500, rank = TRUE, conf = 0.95)
coeff_pcc_Th17g <- pcc(total_data_Th17g[,1:38], total_data_Th17g[,39], nboot = 500, rank = TRUE, conf = 0.95)
coeff_pcc_iTregg <- pcc(total_data_iTregg[,1:38], total_data_iTregg[,39], nboot = 500, rank = TRUE, conf = 0.95)


write.csv(coeff_pcc_HPLP$PRCC, "HPLP_prcc.csv")
write.csv(coeff_pcc_Th1$PRCC, "Th1_prcc.csv")
write.csv(coeff_pcc_Th17$PRCC, "Th17_prcc.csv")
write.csv(coeff_pcc_iTreg$PRCC, "iTreg_prcc.csv")
write.csv(coeff_pcc_Tr$PRCC, "Tr_prcc.csv")
write.csv(coeff_pcc_Mreg$PRCC, "Mreg_prcc.csv")
write.csv(coeff_pcc_Mres$PRCC, "Mres_prcc.csv")
write.csv(coeff_pcc_Minf$PRCC, "Minf_prcc.csv")##Plotting PCC graph ###
write.csv(coeff_pcc_eDC$PRCC, "eDC_prcc.csv")
write.csv(coeff_pcc_tDC$PRCC, "tDC_prcc.csv")
write.csv(coeff_pcc_Th1g$PRCC, "Th1g_prcc.csv")
write.csv(coeff_pcc_Th17g$PRCC, "Th17g_prcc.csv")
write.csv(coeff_pcc_iTregg$PRCC, "iTregg_prcc.csv")


##1.
label <- c("epiinfbctdam", "epiTh1dam", "epiTh17dam", "Epiprolifer.","Epicelldeath","EpiIL10h","nTrep","nTdeath","allTrep","iTregtoTh17","Th17toiTreg","nTtoTr","nTtoiTreg", "nTtoTh17","Th1death","Th17death","iTregdeath","Trdeath", "IL10Tr","dummy","Bactkill","Bactdeath", "HPdeathduetoTcells","HPdeath","DCdeath","Monocytedeath","Resmacdeath","TrmacKill","MregDiff","ResmMacrep","Monorep","IFNg","IL10","IL17","IL21","IL6","TGFb","IL12")
x11()
par(las=2)
barplot(coeff_pcc_HPLP$PRCC$original, horiz=TRUE, col=topo.colors(12), xlim = c(-0.6, 0.6), names.arg = label, mgp=c(0.001,0.00001,.0001), las=1)
##2. 
barplot(coeff_pcc_Mres$PRCC$original, horiz=TRUE, col=topo.colors(12), xlim = c(-0.6, 0.6), names.arg = label, mgp=c(0.001,0.00001,.0001), las=1)
barplot(coeff_pcc_Mreg$PRCC$original, horiz=TRUE, col=topo.colors(12), xlim = c(-0.6, 0.6), names.arg = label, mgp=c(0.001,0.00001,.0001), las=1)
barplot(coeff_pcc_Minf$PRCC$original, horiz=TRUE, col=topo.colors(12), xlim = c(-0.6, 0.6), names.arg = label, mgp=c(0.001,0.00001,.0001), las=1)
barplot(coeff_pcc_eDC$PRCC$original, horiz=TRUE, col=topo.colors(12), xlim = c(-0.6, 0.6), names.arg = label, mgp=c(0.001,0.00001,.0001), las=1)
barplot(coeff_pcc_tDC$PRCC$original, horiz=TRUE, col=topo.colors(12), xlim = c(-0.6, 0.6), names.arg = label, mgp=c(0.001,0.00001,.0001), las=1)
barplot(coeff_pcc_Tr$PRCC$original, horiz=TRUE, col=topo.colors(12), xlim = c(-0.6, 0.6), names.arg = label, mgp=c(0.001,0.00001,.0001), las=1)
barplot(coeff_pcc_Th1$PRCC$original, horiz=TRUE, col=topo.colors(12), xlim = c(-0.6, 0.6), names.arg = label, mgp=c(0.001,0.00001,.0001), las=1)
barplot(coeff_pcc_Th17$PRCC$original, horiz=TRUE, col=topo.colors(12), xlim = c(-0.6, 0.6), names.arg = label, mgp=c(0.001,0.00001,.0001), las=1)
barplot(coeff_pcc_iTreg$PRCC$original, horiz=TRUE, col=topo.colors(12), xlim = c(-0.6, 0.6), names.arg = label, mgp=c(0.001,0.00001,.0001), las=1)
barplot(coeff_pcc_Th1g$PRCC$original, horiz=TRUE, col=topo.colors(12), xlim = c(-0.6, 0.6), names.arg = label, mgp=c(0.001,0.00001,.0001), las=1)
barplot(coeff_pcc_Th17g$PRCC$original, horiz=TRUE, col=topo.colors(12), xlim = c(-0.6, 0.6), names.arg = label, mgp=c(0.001,0.00001,.0001), las=1)
barplot(coeff_pcc_iTregg$PRCC$original, horiz=TRUE, col=topo.colors(12), xlim = c(-0.6, 0.6), names.arg = label, mgp=c(0.001,0.00001,.0001), las=1)


###### PLotting PRCC ######
summary <- print(coeff_pcc_Minf)
plot(summary$original, ylim=c(-1,1), xlab='', ylab='Coefficient',axes=FALSE)
axis(2)

axis(1, at=seq(1:38), labels=label, las=2, cex = 1)
par(ps = 8, cex = 1, cex.main = 1)
mtext(text='Parameter', side=1, line=6)
box()
for(i in 1:38) lines(c(i,i),c(summary[i,4], summary[i,5]))
abline(h=0, col="blue") # draws a horizontal line

#### P values of PRCC ####
## https://www.rdocumentation.org/packages/epiR/versions/0.9-87/topics/epi.prcc

#install.packages("epiR")

library(epiR)
pvalue_LPHP <- epi.prcc(total_data_HPLP, sided.test = 2)

coeff_prcc_pMreg <- epi.prcc(total_data_Mreg, sided.test = 2)
coeff_prcc_pMres <- epi.prcc(total_data_Mres, sided.test = 2)
coeff_prcc_pMinf <- epi.prcc(total_data_Minf, sided.test = 2)

coeff_prcc_pTh1 <- epi.prcc(total_data_Th1, sided.test = 2)
coeff_prcc_pTh17 <- epi.prcc(total_data_Th17, sided.test = 2)
coeff_prcc_piTreg <- epi.prcc(total_data_iTreg, sided.test = 2)

coeff_prcc_pTr <- epi.prcc(total_data_Tr, sided.test = 2)
coeff_prcc_peDC <- epi.prcc(total_data_eDC, sided.test = 2)
coeff_prcc_ptDC <- epi.prcc(total_data_tDC, sided.test = 2)

coeff_prcc_pTh1g <- epi.prcc(total_data_Th1g, sided.test = 2)
coeff_prcc_pTh17g <- epi.prcc(total_data_Th17g, sided.test = 2)
coeff_prcc_piTregg <- epi.prcc(total_data_iTregg, sided.test = 2)


write.csv(pvalue_LPHP, "HPLP_p.csv")
write.csv(coeff_prcc_pTh1, "Th1_p.csv")
write.csv(coeff_prcc_pTh17, "Th17_p.csv")
write.csv(coeff_prcc_piTreg, "iTreg_p.csv")
write.csv(coeff_prcc_pTr, "Tr_p.csv")
write.csv(coeff_prcc_pMreg, "Mreg_p.csv")
write.csv(coeff_prcc_pMres, "Mres_p.csv")
write.csv(coeff_prcc_pMinf, "Minf_p.csv")##Plotting PCC graph ###
write.csv(coeff_prcc_peDC, "eDC_p.csv")
write.csv(coeff_prcc_ptDC, "tDC_p.csv")
write.csv(coeff_prcc_pTh1g, "Th1g_p.csv")
write.csv(coeff_prcc_pTh17g, "Th17g_p.csv")
write.csv(coeff_prcc_piTregg, "iTregg_p.csv")

######## Pvalues #########

HPLP <- which(pvalue_LPHP$p.value < 0.05) # 8 11 12 25 29 34
Mreg <- which(coeff_prcc_pMreg$p.value < 0.05) #4, 12, 14 23, 27, 35
Mres <- which(coeff_prcc_pMres$p.value < 0.05) # 16, 24, 26, 29
Minf <- which(coeff_prcc_pMinf$p.value < 0.05) # 1, 4, 6, 8, 23, 25, 27, 35
Th1 <- which(coeff_prcc_pTh1$p.value < 0.05) #empty
Th17 <- which(coeff_prcc_pTh17$p.value < 0.05) #empty
iTreg <- which(coeff_prcc_piTreg$p.value < 0.05) #empty
Tr <- which(coeff_prcc_pTr$p.value < 0.05)
eDC <- which(coeff_prcc_peDC$p.value < 0.05)
tDC <- which(coeff_prcc_ptDC$p.value < 0.05)
Th1g <- which(coeff_prcc_pTh1g$p.value < 0.05) #empty
Th17g <- which(coeff_prcc_pTh17g$p.value < 0.05) #empty
iTregg <- which(coeff_prcc_piTregg$p.value < 0.05) #empty

write.csv(HPLP, "HPLP.csv")
write.csv(Th1, "Th1.csv")
write.csv(Th17, "Th17.csv")
write.csv(iTreg, "iTreg.csv")
write.csv(Tr, "Tr.csv")
write.csv(Mreg, "Mreg.csv")
write.csv(Mres, "Mres.csv")
write.csv(Minf, "Minf.csv")
write.csv(eDC, "eDC.csv")
write.csv(tDC, "tDC.csv")
write.csv(Th1g, "Th1g.csv")
write.csv(Th17g, "Th17g.csv")
write.csv(iTregg, "iTregg.csv")

##Plotting PCC graph ###

##### Scatter plots ######

## LP HP ###
library("multiplot")
library("Rmisc")
library(ggplot2)

# p8 <- ggplot(total_data_HPLP, aes(x=total_data_HPLP$X8, y=total_data_HPLP$X38)) +
#   geom_point(alpha = 0.4, colour = "black")+
#   geom_smooth(method=lm, color = 'black', alpha = 0.4)
# p8_title <- paste("[PRCC, p-value] =", coeff_pcc_HPLP$PRCC$original[8], pvalue_LPHP$p.value[8])
# p8 <- p8+ggtitle(paste0(p8_title,"\n"))
# 
# p11 <- ggplot(total_data_HPLP, aes(x=total_data_HPLP$X11, y=total_data_HPLP$X38)) +
#   geom_point(alpha = 0.4, colour = "black")+
#   geom_smooth(method=lm, color = 'black', alpha = 0.4)
# p11_title <- paste("[PRCC, p-value] =", coeff_pcc_HPLP$PRCC$original[9], pvalue_LPHP$p.value[11])
# p11 <- p11+ggtitle(paste0(p11_title,"\n"))
# 
# p12 <- ggplot(total_data_HPLP, aes(x=total_data_HPLP$X12, y=total_data_HPLP$X38)) +
#   geom_point(alpha = 0.4, colour = "black")+
#   geom_smooth(method=lm, color = 'black', alpha = 0.4)
# p12_title <- paste("[PRCC, p-value] =", coeff_pcc_HPLP$PRCC$original[12], pvalue_LPHP$p.value[12])
# p12<- p12+ggtitle(paste0(p12_title,"\n"))
# 
# p25 <- ggplot(total_data_HPLP, aes(x=total_data_HPLP$X25, y=total_data_HPLP$X38)) +
#   geom_point(alpha = 0.4, colour = "black")+
#   geom_smooth(method=lm, color = 'black', alpha = 0.4)
# p25_title <- paste("[PRCC, p-value] =", coeff_pcc_HPLP$PRCC$original[25], pvalue_LPHP$p.value[25])
# p25<- p25+ggtitle(paste0(p25_title,"\n"))
# 
# p29 <- ggplot(total_data_HPLP, aes(x=total_data_HPLP$X29, y=total_data_HPLP$X38)) +
#   geom_point(alpha = 0.4, colour = "black")+
#   geom_smooth(method=lm, color = 'black', alpha = 0.4)
# p29_title <- paste("[PRCC, p-value] =", coeff_pcc_HPLP$PRCC$original[29], pvalue_LPHP$p.value[29])
# p29<- p29+ggtitle(paste0(p29_title,"\n"))
# 
# 
# p34 <- ggplot(total_data_HPLP, aes(x=total_data_HPLP$X34, y=total_data_HPLP$X38)) +
#   geom_point(alpha = 0.4, colour = "black")+
#   geom_smooth(method=lm, color = 'black', alpha = 0.4)
# p34_title <- paste("[PRCC, p-value] =", coeff_pcc_HPLP$PRCC$original[34], pvalue_LPHP$p.value[34])
# p34<- p34+ggtitle(paste0(p34_title,"\n"))
# 
# #multiplot(p8, p11, p12, p25, p29, p34, cols = 2)
# library("gridExtra")
# grid.arrange(p8,p11,p12,p25, p29, p34, ncol = 2, main = "H. pylori in LP- parameters")
# library(easyGgplot2)
# ggplot2.multiplot(p8,p11,p12,p25, p29, p34, col = 2)


### Mresident ####
p16 <- ggplot(total_data_Mres, aes(x=total_data_Mres$X16, y=total_data_Mres$X38)) +
  geom_point(alpha = 0.4, colour = "black")+
  geom_smooth(method=lm, color = 'black', alpha = 0.4)
p16_title <- paste("[PRCC, p-value] =", coeff_pcc_Mres$PRCC$original[16], coeff_prcc_pMres$p.value[8])
p16 <- p16+ggtitle(paste0(p16_title,"\n"))

library(easyGgplot2)
ggplot2.multiplot(p16,p24,p26, p29, col = 2)

