####
#Introduction 

#I am interested in the comparisons of how trees can differ in terms of congruence between different genes of the same family, Elapidae. Both COI and Cytochrome b are mitochondrial genes and it is interesting to see how different genes with similar functions can affect data of similar species. COI (Cytochrome c oxidase subunit I) is a gene that is involved in molecular variance within populations but lower among populations, however Cytochrome b is the opposite, which will be interesting to see how this affects the data (Liu et al., 2020). It is known that phylogenetic trees are used for evolutionary biology, and recently a congruence index has been created to help identify topological similarity in these trees (de Vienne et al., 2007). I would like to see how topologically similar two genetic trees are for different genes and why. Further, I want to delve into model building and discover what models would be the best fit when building these trees, using tools such as maximum likelihood estimation. Maximum likelihood estimation are useful for finding the most accurate trees in terms of evolution based on the data that is being used and will be a very interesting tool to build my understanding of models (Yuan et al., 2013). 

#In terms of sub-area of research or gaps in knowledge a growing study method is measuring congruence using R tools. Majority of the methods that currently exist are designed to measure heterogeneity and incongruence for data that is not typically incongruent (Leigh, 2011). One newer tool that has been used in recent studies is CADM, Congruence Among Distance Matrices, a test that estimates the level and can improve type I errors (Campbell, 2011). This is an incredibly useful tool; however it is not useful in my specific study due to how the format must be. This can create gaps in knowledge when the tools can not be widely applied to all datasets and I believe this study will be a great way of comparing models using existing tools and explaining why the need for a proper tool is imminent. Another sub-area of research is co-phylogeny, a visual way of measuring differences in models based on the trees I can construct. This plot can put both trees side by side with lines connecting species, allowing me to notice patterns in clustering, and in turn, evolution.

#The specific objective of my study will be to find the best model, and in turn tree, in R that can measure congruence of both COI and cytochrome b in Elapidae, a family of snakes. In order to do this I will expand on the topics and code I used for assignment 2 in order to improve my skills and better comprehend why modelling decisions are made and how to pick the best one for my data. Using model comparative techniques, I will find the model of best fit and then compare the phylogenetic trees for multiple alignments. The goal of this project is to use these models of best fit to then create a topological congruence plot for both COI and Cytochrome b. This study will answer questions about co-phylogeny but will mainly be about finding the best tools to create the most accurate trees, including the congruence and phylogeny trees.

####
#Description of Data set 

#For this data set I decided to download all of the data available for Elapidae and the COI gene and then all of the data available for Elapidae and Cytochrome b. This searched through the database NCBI to find the information needed and then combined it into one data set of sequences which I could then use for alignment. I first obtained the data on December 1, but in my final check of the dataset I redownloaded everything on December 16, 2022 in order to ensure it was as up to date as possible with every sequence available. I used the entrez_search function and the database db=“nuccore” in order to retrieve it. The data set was converted intro a stringset and then selected to include only the title and the sequence. The title was very long and included the species name within it so I isolated for that and created a new column for identification purposes. For both of the COI gene and Cyt b gene there were initially 150 samples but I chose to go further with only unique species. In terms of their nature, Cytochrome b and COI are the most commonly used genes for fields of study, including forensic studies, and due to their common use for sequencing I decided they would be the best in terms of comparison studies as they would be readily available (Nevo, 2001). The data set itself is not from one study and can not be cited, however the website that contains all of the data is NCBI. 

####
#Code Section 1. Data Acquisition, Exploration, Filtering, and Quality Control

#load all of the packages that will be needed for the script to keep them in one place
library(rentrez)
library(BiocManager)
BiocManager::install("Biostrings")
#no updates
n
library(Biostrings)
library(stringr)
library(tidyverse)
library(DECIPHER)
library(dendextend)
library(ape)
library(phangorn)
library(phytools)

#Define the amount of missing data to accept for this analysis
missing.data <- 0.01

#Define the amount of sequence length variability to accept
length.var <- 50

#search for data from the Nucleotide database, I am looking specifically for snakes (cobras, coral snakes, sea snakes, and a few other kinds) and the cytochrome b gene.
load.cytb<-entrez_search(db = "nuccore", term="(Elapidae[ORGN] AND CytB[Gene] AND 100:1500[SLEN])", retmax = 150) 
summary(load.cytb)

#ensure that there are 150 
length(load.cytb$ids)  

#repeat the search for results with only the COI gene
load.COI<-entrez_search(db = "nuccore", term="(Elapidae[ORGN] AND COI[Gene] AND 100:1500[SLEN])", retmax = 150) 
summary(load.COI)

#ensure again that this worked
length(load.COI$ids)  

#the counts should be the same as when you click the variable for the two load variables, which they are
load.cytb$count
load.COI$count

#load the cytochrome b nucleotide sequences in FASTA format.
snake_cytb <- entrez_fetch(db = "nuccore", id = load.cytb$ids, rettype = "fasta")

#save this onto my own computer for reference with no changes to the original file so I can keep a copy of the original dataset
write(snake_cytb, "snake_cytb.fasta", sep = "\n")

#repeat for COI 
snake_COI <- entrez_fetch(db = "nuccore", id = load.COI$ids, rettype = "fasta")
write(snake_COI, "snake_COI.fasta", sep = "\n")

#save these into our working environment to continue the data manipulation phase
cytb.string<-readDNAStringSet("snake_cytb.fasta")
COI.string<-readDNAStringSet("snake_COI.fasta")

#convert the format into a dataframe to make it workable
data.cytb<-data.frame(Title = names(cytb.string), Sequence = paste(cytb.string))
data.COI<-data.frame(Title = names(COI.string), Sequence = paste(COI.string))

#add a column that will include the Species of each sequence. This is done to help visualize and understand the data
data.cytb$Species<-word(data.cytb$Title, 2L, 3L)
data.cytb<-data.cytb[, c("Title", "Species", "Sequence")]

data.COI$Species<-word(data.COI$Title, 2L, 3L)
data.COI<-data.COI[, c("Title", "Species", "Sequence")]

#my next step will be to filter to a manageable amount of species, I will start by checking how many unique species I have for each search.

length(unique(data.cytb$Species))
length(unique(data.COI$Species))

#There are 36 unique species for the cytochrome b and 30 for COI. 

#ensure the results always come out the same
set.seed(190)

#I will filter the data now, only using unique species, removing the NAs in the process along with the gaps. I also set it to be within 50 base length of each other and 1% N's

subs_cytb <- data.cytb %>%
  group_by(Species) %>% 
  filter(!is.na(Sequence)) %>%
  mutate(Nucleotides = str_remove_all(Sequence, "^N+|N+$|-")) %>%
  filter(str_count(Nucleotides, "N") <= (missing.data * str_count(Nucleotides))) %>%
  filter(str_length(Nucleotides) >= median(str_length(Nucleotides)) - length.var & str_length(Nucleotides) <= median(str_length(Nucleotides)) + length.var) %>%
  sample_n(1)

subs_COI <- data.COI %>% 
  group_by(Species) %>% 
  filter(!is.na(Sequence)) %>%
  mutate(Nucleotides = str_remove_all(Sequence, "^N+|N+$|-")) %>%
  filter(str_count(Nucleotides, "N") <= (missing.data * str_count(Nucleotides))) %>%
  filter(str_length(Nucleotides) >= median(str_length(Nucleotides)) - length.var & str_length(Nucleotides) <= median(str_length(Nucleotides)) + length.var) %>%
  sample_n(1)

#This narrowed it down to 36 and 30 species respectively, a manageable amount for data comparisons

#a quick failsafe to make sure this worked, the subset I just made should be equal to the amount of unique species from the variable I used to make these subsets.
all.equal(length(unique(data.cytb$Species)), nrow(subs_cytb))
all.equal(length(unique(data.COI$Species)), nrow(subs_COI))

#add a column counting the length for making a boxplot
subs_cytb$length<-str_count(subs_cytb$Nucleotides)
subs_COI$length<-str_count(subs_COI$Nucleotides)

#allow for two plots to be side by side for comparison
par(mfrow=c(1,2),mar = c(4, 4, 4, 4))

#visualize a boxplot for the length of the nucelotide sequences
boxplot(subs_COI$length,ylab="Length of Sequence",main="COI",col="light blue")
boxplot(subs_cytb$length,ylab="Length of Sequence",main="Cyt b",col="light yellow")
mtext("Figure 1. Boxplots of both Cytochrome b and COI genes for sequence length to identify any outliers and compare medians",side=1,outer=TRUE,cex=1)

#there are quite a few outliers on both graphs, telling me there could still be some more sequences skewing the data/alignments

#let's take a look at the medians and compare for both genes
median(subs_cytb$length)
median(subs_COI$length)

#Some summary code to use to better understand the data I am working with. First, check the unique species I just selected for 
unique(subs_cytb$Species)
unique(subs_COI$Species)

#ensure the sum of NAs are 0
sum(is.na(subs_cytb$Nucleotides))
sum(is.na(subs_COI$Nucleotides))

#check to make sure the gaps have been removed, should both be 0
sum(str_count(subs_cytb$Nucleotides, "-"))
sum(str_count(subs_COI$Nucleotides, "-"))

#use the summary function to see what I am working with in terms of median, mean, and quantiles 
summary(str_count(subs_cytb$Nucleotides))
summary(str_count(subs_COI$Nucleotides))

#Reformat into "DNAStringSet". First make it a data frame and then convert it to the proper format for alignment. 
cytb_sub<-as.data.frame(subs_cytb)
COI_sub<-as.data.frame(subs_COI)
cytb_sub$Nucleotides2<-DNAStringSet(cytb_sub$Nucleotides)
COI_sub$Nucleotides2<-DNAStringSet(COI_sub$Nucleotides)

#Ensure it now says "DNAStringSet"
class(cytb_sub$Nucleotides2)
class(COI_sub$Nucleotides2)

#I will now take a look at the sequences, see if I can find any patterns before doing any alignment
BrowseSeqs(cytb_sub$Nucleotides2)
BrowseSeqs(COI_sub$Nucleotides2)

#I can see that there are a few shorter sequences in the cytb nucleotides, I can also see a large chunk of the COI seems to be somewhat nicely aligned already. I will keep these in mind as I made my first few visualization plots, now I have an idea of what I am looking for.

#Layer the histograms on top of each other to visualize the two genes side by side
par(mfrow=c(1,1),oma = c(2, 0, 0, 0))
hist(nchar(COI_sub$Nucleotides2),xlab = "Nucleotide Length",col=rgb(1,0,0,0.2),xlim=c(300,1000))
hist(nchar(cytb_sub$Nucleotides2),xlab = "Nucleotide Length",main="",col=rgb(0,0,1,0.2),xlim=c(300,1000),add=TRUE)
legend('topright', c('Cytochrome b', 'COI','Both'),
       fill=c(rgb(0,0,1,0.2), rgb(1,0,0,0.2),rgb(0.6,0,0.4,0.3)))
mtext("Figure 2. Histogram of Cytochrome b and COI genes for length to get a better look at the distribution",side=1,outer=TRUE,cex=1)

#The line of code did not correctly remove my sequences within a distance of the median, I will manually do this for the second alignment and discuss this is in my writeup. 

#I will align the sequences prior to deciding which sequences to remove, as it is not completely clear yet how it is affecting the data

#Create the identifier using the species names
names(cytb_sub$Nucleotides2)<-cytb_sub$Species
names(COI_sub$Nucleotides2)<-COI_sub$Species

#Using the default settings, try an alignment to see what (if anything) needs to change
align1.cytb <- DNAStringSet(muscle::muscle(cytb_sub$Nucleotides2, log = "log.tx", verbose = T), use.names = T)
align1.COI <- DNAStringSet(muscle::muscle(COI_sub$Nucleotides2, log = "log.tx", verbose = T), use.names = T)

#check how they look
BrowseSeqs(align1.cytb)
BrowseSeqs(align1.COI)

#for cytochrome b there are two sequences with under 350 nucleotides, I will most likely remove those and try for a better alignment.There are a lot of significant gaps in those two sequences right in the middle of the alignment.
#for COI I am going to most likely remove Ophiophagus hannah, this species appears to be the biggest contributor to changing the alignment and can be removed without making the study less effective.

#the best method of ensuring I am removing the correct sequences is to create a dendrogram. 

#find the mean of all the gap counts and get a better understanding of the alignment
summary(unlist(lapply(align1.cytb, str_count, "-")))
summary(unlist(lapply(align1.COI, str_count, "-")))

#the mean of the gaps is quite a lot higher in the cytochrome b alignment versus the COI alignment.

#For the last element of visualizing the data, I will create a dendrogram by clustering the data and creating the tree.

#format as DNAbin to use for the distance matrix
dnaBIN.cytb1<-as.DNAbin(align1.cytb)
dnaBIN.COI1<-as.DNAbin(align1.COI)
class(dnaBIN.cytb1)
class(dnaBIN.COI1)

#choose model and clustering method. Since I am using COI, I will choose the JC method (Evolution and Genomics, 2022).This is a starting point, I will do further tests on models in the main code section.
chosen.model<-"JC69"

#I tried many different methods, single clearly did not work, and ward.D is great for continous data and is why I made this decision
clustering.method<-"ward.D"

#build a distance matrix using the chosen model
distMat.cytb<-dist.dna(dnaBIN.cytb1, model=chosen.model,as.matrix = TRUE, pairwise.deletion = TRUE)
distMat.COI<-dist.dna(dnaBIN.COI1, model=chosen.model,as.matrix = TRUE, pairwise.deletion = TRUE)

#cluster it according to the clustering method
distMat.cytb=as.dist(distMat.cytb)

#Code incorporated from: https://stackoverflow.com/questions/54153670/remove-na-values-from-distance-matrix-in-r  (cited below)
#remove Nas and then remake the matrix
x = as.matrix(distMat.cytb)
x = x[rowSums(is.na(x)) == 0, colSums(is.na(x)) == 0, drop = FALSE]
p=as.dist(x)

#cluster based off the predetermined method
cytb.cluster<-hclust(p,method=clustering.method)

distMat.COI=as.dist(distMat.COI)
COI.cluster<-hclust(distMat.COI,method=clustering.method)

#Now that the code has been checked, I will remove the sequences from cytochrome b that are below 500 and above 900 to make the lengths closer together
subs2_cytb<-subs_cytb[-c(1:7,17,20:22,25,34:35),]

#For COI I will remove the sequences below and including a length of 550
subs2_COI<-subs_COI[-c(5,26,28:29),]

#change format again, using the same steps as above
cytb_sub2<-as.data.frame(subs2_cytb)
COI_sub2<-as.data.frame(subs2_COI)
cytb_sub2$Nucleotides2<-DNAStringSet(cytb_sub2$Nucleotides)
COI_sub2$Nucleotides2<-DNAStringSet(COI_sub2$Nucleotides)

names(cytb_sub2$Nucleotides2)<-cytb_sub2$Species
names(COI_sub2$Nucleotides2)<-COI_sub2$Species

#Repeat alignments for these updated datasets
align2.cytb <- DNAStringSet(muscle::muscle(cytb_sub2$Nucleotides2, log = "log.tx", verbose = T), use.names = T)
align2.COI <- DNAStringSet(muscle::muscle(COI_sub2$Nucleotides2, log = "log.tx", verbose = T), use.names = T)

#view the new alignment
BrowseSeqs(align2.cytb)
BrowseSeqs(align2.COI)

#I can clearly see that for cytb, the third sequence is negatively affecting the alignment and I will remove that as well. COI alignment 2 looks good and ready for further analysis.

subs3_cytb<-subs2_cytb[-3,]
cytb_sub3<-as.data.frame(subs3_cytb)
cytb_sub3$Nucleotides2<-DNAStringSet(cytb_sub3$Nucleotides)
names(cytb_sub3$Nucleotides2)<-cytb_sub3$Species
align3.cytb <- DNAStringSet(muscle::muscle(cytb_sub3$Nucleotides2, log = "log.tx", verbose = T), use.names = T)
BrowseSeqs(align3.cytb)

#This alignment looks much better. Repeat the steps and compare dendrograms

dnaBIN.cytb2<-as.DNAbin(align3.cytb)
dnaBIN.COI2<-as.DNAbin(align2.COI)

distMat.cytb2<-dist.dna(dnaBIN.cytb2, model=chosen.model,as.matrix = TRUE, pairwise.deletion = TRUE)
distMat.COI2<-dist.dna(dnaBIN.COI2, model=chosen.model,as.matrix = TRUE, pairwise.deletion = TRUE)

#cluster it according to the clustering method
par(mfrow=c(2,2),mar = c(3, 3, 3, 3))
distMat.cytb2=as.dist(distMat.cytb2)

#There are no Na's this time! No need to alter the distance matrix

#cluster based off the predetermined method
cytb.cluster2<-hclust(distMat.cytb2,method=clustering.method)

#create the dendrogram and colour code the clusters for the original alignments
dend<-as.dendrogram(cytb.cluster)
de=color_branches(dend,k=6)%>%
  set("labels_cex",0.8)%>%
  set("branches_lwd", 2)
plot(de,mar=c(.1,.1,2,.1),ylim=c(-0.1,0.25),main="Cyt b")

dendr<-as.dendrogram(COI.cluster)
de2=color_branches(dendr,k=8)%>%
  set("labels_cex",0.8)%>%
  set("branches_lwd", 2)
plot(de2,mar=c(.1,.1,2,.1),ylim=c(-0.1,0.25),main="COI")

#create the dendrogram and colour code the clusters to visualize the final alignments
dend.cytb<-as.dendrogram(cytb.cluster2)
de2=color_branches(dend.cytb,k=6)%>%
  set("labels_cex",0.8)%>%
  set("branches_lwd", 2)
plot(de2,mar=c(.1,.1,2,.1),ylim=c(-0.1,0.25),main="Cyt b - Final")

distMat.COI2=as.dist(distMat.COI2)
COI.cluster2<-hclust(distMat.COI2,method=clustering.method)
dendr.COI<-as.dendrogram(COI.cluster2)
de.COI=color_branches(dendr.COI,k=8)%>%
  set("labels_cex",0.8)%>%
  set("branches_lwd", 2)
plot(de.COI,mar=c(.1,.1,2,.1),ylim=c(-0.1,0.25),main="COI - Final")
mtext("Figure 3. Cluster Dendrograms (hierarchical clustering) of both Cytochrome b and COI genes from alignment 1",side=1,outer=TRUE,cex=0.9)

####
#Main Software Tools

#I filtered and manipulated the data described above to prepare it for alignment and tried out a few alignments to find the one of best fit. Next, I will use the matrices to try some maximum likelihood estimates and pick a model that best fits the data to create the most accurate tree. Some examples include the SH (Shimodaira-Hasegawa) test and the Kruskal-Wallis Test. One of my more novel plots will be a plot of both multidimensional scaling (MDS) for both models (Cyt b and COI), to compare the similarity and differences for both models and see how they fit comparatively. This decision was made to expand my knowledge on congruence and also what is actually occurring in the models and why I chose them. I considered many alternatives, the biggest of which is CADM.global, a method of measuring congruence. Unfortunately, I could not get the function to correctly compare my models due to the specific .txt format it required. I calculated bootstraps and plotted this onto my trees of the final models, as a way of seeing how replicable these trees were. Finally, I will create a topological congruence graph after performing all those steps. I built off of the vignette for phangorn in order to understand maximum likelihood, but I added my own plots and model comparison methods in order to improve this for my specific study.

####
#Code section 2. Main Analysis

#convert the alignment into a matrix
align.cyt<-as.matrix(align3.cytb)
align.COI<-as.matrix(align2.COI)

#convert the matrix into phyDat format for maximum likelihood analysis
formt.cyt<-as.phyDat(align.cyt)
formt.COI<-as.phyDat(align.COI)

#make a distance matrix
dm1<-dist.ml(formt.cyt)
dm2<-dist.ml(formt.COI)

#tree 1/2 is the same tree as above, just clustered in a different order 
tree1<-upgma(dm1)
tree2<-upgma(dm2)

#run a model test using common models to check which one would be the best
mt <- modelTest(formt.cyt, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR"), control = pml.control(trace = 0))
mt
mt2 <- modelTest(formt.COI, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR"), control = pml.control(trace = 0))
mt2

#BIC is useful for model selection based on existing data, will choose a model based on highest BIC:
fit_cyt <- pml_bb(mt, control = pml.control(trace = 0))
fit_cyt
fit_COI <- pml_bb(mt2, control = pml.control(trace = 0))
fit_COI

#Let's compare this to the original model we chose above, JC.
fit.cyt<-pml(tree1,data=formt.cyt)
fit.COI<-pml(tree2,data=formt.COI)

#Use this basic model fit and update it for the JC method. Use NNI for trees that are one nearest neighbor interchange away
JC.cyt<-optim.pml(fit.cyt, rearrangement="NNI")
JC.COI<-optim.pml(fit.COI, rearrangement="NNI")

#Compare these two models:

#Use BIC method
BIC(fit_cyt)
BIC(JC.cyt)

BIC(fit_COI)
BIC(JC.COI)

#Calculating formal metrics of phylogenetic congruence in model determination:

#Shimodaira-Hasegawa test: can compare the congruence of two trees.
SH.test(JC.cyt,fit_cyt)
SH.test(JC.COI,fit_COI)
#P-value is large, big difference in log-likelihood, we can assume these trees are congruent

#going beyond this, a Kruskal test will tell me the Kruskal stress statistic which can tell me the discrepancy between the matrices. A low statistic = higher congruence.

#first, put the trees into distance matrix from both models
JC.mat<-cophenetic(JC.cyt$tree)
fit.mat<-cophenetic(fit_cyt$tree)
JC.mat.COI<-cophenetic(JC.COI$tree)

congr<-kruskal.test(JC.mat,fit.mat)   #this also gives a p-value which is very small
congr$statistic

#the small p-value indicates the two distance matrices are significantly correlated

#Now will try out MDS, multidimensional scaling on both matrices and the plot them together to visually look for differences.
#plot out the MDS for both distance matrices
par(mfrow=c(1,1),mar = c(4, 4, 4, 4))
multi1<-cmdscale(JC.mat)
multi3<-cmdscale(JC.mat.COI)
both<-rbind(multi1,multi3)
plot(both,col=c("blue","black"),pch=15)
legend('topleft', c('Cytochrome b', 'COI'),
       fill=c("blue","black"))
mtext("Figure 4. Multidimensional scaling for JC models for both Cytochrome b and COI",side=1,outer=TRUE,cex=1)

#The JC models appear to be a better fit, though the difference is not large, and we will move forward with these.
#let's compare these two models to each other - the two different genes
anova(JC.cyt,JC.COI)


#JC.cyt is a better fit for its data than JC.COI is.

#Calculate how well the edges of the tree are supported, this is done using bootstrap
bs.cyt<-bootstrap.pml(JC.cyt,bs=100, optNni=TRUE, control = pml.control(trace = 0))
bs.COI<-bootstrap.pml(JC.COI,bs=100, optNni=TRUE, control = pml.control(trace = 0))

par(mfrow=c(1,2),mar = c(2, 2, 2, 2))

plotBS(midpoint(JC.cyt$tree),bs.cyt,p=50,type="p",bs.col="blue",main="Cyt b")
plotBS(midpoint(JC.COI$tree),bs.COI,p=50,type="p",bs.col="red",main="COI")
mtext("Figure 5. Plot using the JC model and bootsrap support values",side=1,outer=TRUE,cex=1)


#create the Co-phylogenetic plot and add lines showing the same species
par(mfrow=c(1,1),mar = c(0.5, 0.5, 0.5, 0.5))
comb<-cophylo(JC.cyt$tree,JC.COI$tree,print=TRUE)
plot(comb,mar=c(.1,.1,2,.1),scale.bar=c(.005,.05),ylim=c(-.2,1),link.type="curved",link.lwd=4,link.lty="solid",link.col=make.transparent("red",0.25),fsize=0.9)
title("Cytochrome b                                                        COI")
mtext("Figure 6. Topological congruence on a co-phylogenetic plot",side=1,outer=TRUE,cex=1)

####
#Results and Discussion

#My original plan was to measure congruence for COI and Cytochrome b, however this was more complicated than I anticipated. I was able to create a plot, figure 4, that plotted both congruences and I could clearly see some differences. There were some similar patterns in terms of MDS distributions, but there were some more extreme plot points for COI than there were for Cytochrome b. I also was able to answer my main question on how these two genes differ in terms of co-phylogeny and evolution in figure 6, where I noticed some patterns in terms of cluster changing. Naja kaouthia and Naja melanoleuca are clustered together for Cytochrome b but have quite a few clusters in between them for COI. It would be really interesting as a future direction to expand this project with much more species, being able to see more variety in terms of clustering between the two models when similar species are used. The results were somewhat expected, I did anticipate seeing cluster changes, however in terms of congruence I anticipated being able to compare the models more accurately, which also would have been improved if I had a bigger data set of equal size for both genes. Knowing that COI has lower variance among populations, I expected Cytochrome b to be less accurate for this data as it is not from one population (Liu et al., 2020). This is interesting in terms of the MDS plot and the patterns there.

#A key caveat of my study was shown through the plot of the bootstrap values. The bootstrap plot was able to quantify limitations and successes. There were quite a few 100s, which indicate support of the results that when repeating the generation of a phylogenetic tree 100 times we receive the same clustering (Efron et al., 1996). However, these plots indicate that there are likely better models or other factors to consider in the future before settling on these trees. Another key part of this study was that my initial subset did not work when filtering data. I tried to filter only rows that were within 50 of the length of the median nucleotide sequence, but this was not the case. This incorrect filtering could have impacted my data in a way I did not catch, which could introduce inconsistencies to my data. Since I did catch this, I manually removed rows based on my histogram plot which helped with making better alignments. Data availability was not a problem for me as I managed to find 150 sequences for both genes of my specified family, however since I was downloading from NCBI and all their data, there could introduce many kinds of bias. This includes if methodology for sequencing was not standard among different studies and if there were errors. It would be interesting to go back to the data itself and try and sort through it to see what kinds of bias I could find and how I could fix it in a future study.

#The next steps for this research could be to use what we learned about the topological congruence to enhance other studies, such as phylogeographic studies and how geography affects these genes, or evolutionary studies. In terms of my own next steps, with more time or resources I would expand on the congruence section. I think I was on the right track with my model testing techniques, however unfortunately since the two main models for COI and Cyt b had different lengths (amount of species), some of the tests could not be used to compare these models but only to compare models of the same gene. I would like to filter out even more and choose to use the same number of species for both studies and see how these congruence tests would change the results I found. The results finding different patterns in the co-phylogeny plot for some of the clustering is something I would love to follow up on in a greater scale. I would be very interested in seeing how evolution is related to these findings and better understand what they mean. This study was a great tool in terms of developing a relationship between congruence and phylogeny and expanding on this could yield some really interesting scientific discoveries. 

#Something I really wanted to focus on was improving my assignment 2, as I feel like I learned a lot of useful skills and knowledge from that study, but there were some areas of focus I did not feel like I developed enough. In this study I managed to include multiple new tests, including maximum likelihood estimation, and plots, including the MDS plot, in order to dive into the data in a way that taught me much more about the models themselves. This was a great improvement for me instead of just accepting whatever first model I could find. As with the guest speakers, I learned a lot from these hands-on projects in terms of making decisions. I think this was something I showed improvement on between the first couple of assignments and now.


####
#References.

#Campbell, V., Legendre, P., & Lapointe, F.-J. (2011). The performance of the congruence among distance matrices (CADM) test in phylogenetic analysis. BMC Evolutionary Biology, 11(1). https://doi.org/10.1186/1471-2148-11-64 

#de Vienne, D. M., Giraud, T., & Martin, O. C. (2007). A congruence index for testing topological similarity between trees. Bioinformatics, 23(23), 3119–3124. https://doi.org/10.1093/bioinformatics/btm500 

#Efron, B., Halloran, E., & Holmes, S. (1996). Bootstrap confidence levels for phylogenetic trees. Proceedings of the National Academy of Sciences, 93(23), 13429–13429. https://doi.org/10.1073/pnas.93.23.13429 

#Leigh, J. W., Lapointe, F.-J., Lopez, P., & Bapteste, E. (2011). Evaluating phylogenetic congruence in the post-genomic era. Genome Biology and Evolution, 3, 571–587. https://doi.org/10.1093/gbe/evr050 

#Liu, B., Zhang, K., Zhu, K., Shafi, M., Gong, L., Jiang, L., Liu, L., Muhammad, F., & Lü, Z. (2020). Population genetics of konosirus punctatus in Chinese coastal waters inferred from two mtdna genes (Coi and Cytb). Frontiers in Marine Science, 7. https://doi.org/10.3389/fmars.2020.00534 

#Yuan, A., Jeffries, N., Zheng, G., Maloy, S. R., & Hughes, K. T. (2013). Brenner's Encyclopedia of Genetics. Elsevier/Academic Press. p. 328-329

#Nevo, E. Genetic Diversity. In: Encyclopedia of Biodiversity (ed. Levin, S.) p. 16 (Academic Press, 2001).



#Vignettes and Stack Exchange:

#Estimating phylogenetic trees with Phangorn. (2022, September 17). Retrieved December 3, 2022, from https://cran.r-project.org/web/packages/phangorn/vignettes/Trees.html 

#Thomas, G. (2019, February 1). Remove na values from distance matrix in R. Stack Overflow. Retrieved December 5, 2022, from https://stackoverflow.com/questions/54153670/remove-na-values-from-distance-matrix-in-r 
