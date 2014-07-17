#####    #####     #####    #####     #####    #####     #####    #####     #####    #####     #####
# Code to Create the 2011 Area Classification for Output Areas (2011 OAC)
# Feel free to share and reuse with attribution
# Pre-Cluster 1 - Correlation and Within-Cluster Sum of Squares analysis
#####    #####     #####    #####     #####    #####     #####    #####     #####    #####     #####

rm(list = ls()) 
gc()

#OACInputData <- FILE INPUT NAME

##########LIBRARY PACKAGES REQUIRED##########
library(cluster)
library(fpc)
library(maptools)
library(spdep)
library(Hmisc)
library(ggplot2)
library(plotrix)
library(reshape)
library(car)
library(grid)
library(lattice)
library(circular)
library(RColorBrewer)
library(CircStats)
library(hexbin)
library(modeest)
library(plyr)
library(scales)
library(reshape2)
library(MASS)

##########WORKSPACE SELECTION AND DIRECTORY CREATION##########

setwd("/Users/Chris/Documents/2011 OAC")

dir.create("Missing Variables", showWarnings = FALSE)
dir.create("Missing Variables/Variables - All", showWarnings = FALSE)
dir.create("Missing Variables/Variables - Domains", showWarnings = FALSE)
dir.create("Missing Variables/Summary Data", showWarnings = FALSE)
dir.create("Correlation", showWarnings = FALSE)
dir.create("Distribution Plots", showWarnings = FALSE)

##########INPUT##########

#Number of Clusters Range
CN <- 5:9

#Number of k-means loops
KM <- 2

#2011 UK OAC Variable Domain ranges
G1 <- 1:50
G2 <- 51:83
G3 <- 84:103
G4 <- 104:122
G5 <- 123:167

#Variable Domain names
N1 = "Demographic"
N2 = "Household Composition"
N3 = "Housing"
N4 = "Socio-Economic"
N5 = "Employment"

#Variable Domain colours
C1 <- "#E41A1C"
C2 <- "#377EB8"
C3 <- "#4DAF4A"
C4 <- "#984EA3"
C5 <- "#FF7F00"

#Number of columns for plots
FCol <- 9

#Number of rows for plots
FRow <- 7

#Do you wish to produce cluster outputs for the missing variables datasets?  
#Enter "YES"/"NO"
RQMISSOUT<- "NO"

#Do you wish to produce a text file updating you on the progress of the clustering with missing variables datasets?  
#Enter "YES"/"NO"
RQCLUSUP<- "YES"

#Do you have a custom 'VFinal' CSV file you wish to use?
#Enter "YES"/"NO"
RQCUSTOMV <- "NO"

if (RQMISSOUT == "YES")
{
dir.create("Missing Variables/Datasets", showWarnings = FALSE)
}

if (RQCLUSUP == "YES")
{
dir.create("Missing Variables/Cluster Update", showWarnings = FALSE)
}

if (RQCUSTOMV == "YES")
{
vdata<- read.csv(file="VFinal.csv",head=TRUE,sep=",")
}

OACInputData_Input <- read.csv(file="OACInputData.csv",head=TRUE,sep=",")

#OACInputData_Input<- read.csv("D:/2011 UK OAC/Stage 1/CSV Files/OACInputData.csv",head=TRUE,sep=",")


##########INCREASE MEMORY (WINDOWS ONLY##########
if (Sys.info()['sysname']=="Windows")	
{memory.limit(40000)}

##########CORRELATION##########
CorrelationStart <- Sys.time()

OACInputData_Input_CORR <- data.frame(OACInputData_Input, row.names=1)

OACInputData_Input_Variable_Names <-data.frame(colnames(OACInputData_Input_CORR))
OACInputData_Input_Variable_Numbers <-nrow(OACInputData_Input_Variable_Names)
OACInputData_Input_Variable_Number_Range <-data.frame(1:OACInputData_Input_Variable_Numbers)
OACInputData_Input_Variables <- cbind(OACInputData_Input_Variable_Number_Range, OACInputData_Input_Variable_Names)
colnames(OACInputData_Input_Variables) <-c("Variable Number", "Variable Name")

OACInputData_P_CORR <- as.matrix(OACInputData_Input_CORR,rownames.force=TRUE)

OACInputData_P_CORR <- as.matrix(OACInputData_P_CORR)

RowNumber <- nrow(OACInputData_Input_CORR)
ColNumber <- ncol(OACInputData_Input_CORR)

R <- rcorr(OACInputData_P_CORR, type="pearson")$r
p <- rcorr(OACInputData_P_CORR, type="pearson")$P

CORR_STARS <- ifelse(p == 0, "", ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05,
"* ", " "))))

R <- format(round(cbind(rep(-1.11, ncol(OACInputData_P_CORR)), R), 2))[,-1]

OACInputData_P_SIG_CORR <- matrix(paste(R, CORR_STARS, sep=""), ncol=ncol(OACInputData_P_CORR))
diag(OACInputData_P_SIG_CORR) <- paste(diag(R), " ", sep="")
rownames(OACInputData_P_SIG_CORR) <- colnames(OACInputData_P_CORR)
colnames(OACInputData_P_SIG_CORR) <- paste(colnames(OACInputData_P_CORR), "", sep="")

OACInputData_P_CORR_PEARSON_P <- rcorr(OACInputData_P_CORR, type="pearson")$P 
OACInputData_P_CORR_PEARSON_P <- format(round(cbind(rep(-1.11, ncol(OACInputData_P_CORR)), OACInputData_P_CORR_PEARSON_P), 4))[,-1] 
OACInputData_P_SIG_PVALUE <- as.data.frame(OACInputData_P_CORR_PEARSON_P)
OACInputData_P_SIG_CORR <- as.data.frame(OACInputData_P_SIG_CORR)

counta=1
countb=2

OACInputData_P_SEL_PVALUE <-OACInputData_P_SIG_PVALUE[1, 1:ncol(OACInputData_P_SIG_PVALUE)]
OACInputData_P_SEL_FINAL <-OACInputData_P_SEL_PVALUE
names(OACInputData_P_SEL_FINAL)<- c(names(OACInputData_P_SIG_PVALUE))

for(uu in 2:nrow(OACInputData_P_SIG_PVALUE))
	{
	OACInputData_P_SEL_PVALUE <-OACInputData_P_SIG_PVALUE[uu, countb:ncol(OACInputData_P_SIG_PVALUE)]
	OACInputData_P_SEL_CORR <-OACInputData_P_SIG_CORR[uu, 1:counta]

	OACInputData_P_SEL_COMBINE <- cbind(OACInputData_P_SEL_CORR, OACInputData_P_SEL_PVALUE)[1:ColNumber]
	names(OACInputData_P_SEL_COMBINE)<- c(names(OACInputData_P_SIG_PVALUE))
	OACInputData_P_SEL_FINAL <-rbind(OACInputData_P_SEL_FINAL,OACInputData_P_SEL_COMBINE)
	names(OACInputData_P_SEL_FINAL)<- c(names(OACInputData_P_SIG_PVALUE))
	counta=counta+1
	countb=countb+1
	}

OACInputData_S_CORR <- as.matrix(OACInputData_Input_CORR,rownames.force=TRUE)

OACInputData_S_CORR <- as.matrix(OACInputData_S_CORR)

#These calculations may take a few minutes to complete
R <- rcorr(OACInputData_S_CORR, type="spearman")$r
p <- rcorr(OACInputData_S_CORR, type="spearman")$P

CORR_STARS <- ifelse(p == 0, "", ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05,
"* ", " "))))

R <- format(round(cbind(rep(-1.11, ncol(OACInputData_S_CORR)), R), 2))[,-1]

OACInputData_S_SIG_CORR <- matrix(paste(R, CORR_STARS, sep=""), ncol=ncol(OACInputData_S_CORR))
diag(OACInputData_S_SIG_CORR) <- paste(diag(R), " ", sep="")
rownames(OACInputData_S_SIG_CORR) <- colnames(OACInputData_S_CORR)
colnames(OACInputData_S_SIG_CORR) <- paste(colnames(OACInputData_S_CORR), "", sep="")

#These calculations may take a few minutes to complete
OACInputData_S_CORR_SPEARMAN_P <- rcorr(OACInputData_S_CORR, type="spearman")$P 
OACInputData_S_CORR_SPEARMAN_P <- format(round(cbind(rep(-1.11, ncol(OACInputData_S_CORR)), OACInputData_S_CORR_SPEARMAN_P), 4))[,-1] 
OACInputData_S_SIG_PVALUE <- as.data.frame(OACInputData_S_CORR_SPEARMAN_P)
OACInputData_S_SIG_CORR <- as.data.frame(OACInputData_S_SIG_CORR)

counta=1
countb=2

OACInputData_S_SEL_PVALUE <-OACInputData_S_SIG_PVALUE[1, 1:ncol(OACInputData_S_SIG_PVALUE)]
OACInputData_S_SEL_FINAL <-OACInputData_S_SEL_PVALUE
names(OACInputData_S_SEL_FINAL)<- c(names(OACInputData_S_SIG_PVALUE))

for(uu in 2:nrow(OACInputData_S_SIG_PVALUE))
	{
	
	OACInputData_S_SEL_PVALUE <-OACInputData_S_SIG_PVALUE[uu, countb:ncol(OACInputData_S_SIG_PVALUE)]
	OACInputData_S_SEL_CORR <-OACInputData_S_SIG_CORR[uu, 1:counta]

	OACInputData_S_SEL_COMBINE <- cbind(OACInputData_S_SEL_CORR, OACInputData_S_SEL_PVALUE)[1:ColNumber]
	names(OACInputData_S_SEL_COMBINE)<- c(names(OACInputData_S_SIG_PVALUE))
	OACInputData_S_SEL_FINAL <-rbind(OACInputData_S_SEL_FINAL,OACInputData_S_SEL_COMBINE)
	names(OACInputData_S_SEL_FINAL)<- c(names(OACInputData_S_SIG_PVALUE))
	counta=counta+1
	countb=countb+1

	}

OACInputData_CORR <- as.matrix(OACInputData_Input_CORR,rownames.force=TRUE)

OACInputData_CORR <- as.matrix(OACInputData_CORR)

#These calculations may take a few minutes to complete
OACInputData_CORR_PEARSON <- rcorr(OACInputData_CORR, type="pearson") 
OACInputData_CORR_PEARSON_R <- rcorr(OACInputData_CORR, type="pearson")$r 
OACInputData_CORR_PEARSON_P <- rcorr(OACInputData_CORR, type="pearson")$P 
OACInputData_CORR_PEARSON_R <- format(round(cbind(rep(-1.11, ncol(OACInputData_CORR)), OACInputData_CORR_PEARSON_R), 4))[,-1] 
OACInputData_CORR_PEARSON_P <- format(round(cbind(rep(-1.11, ncol(OACInputData_CORR)), OACInputData_CORR_PEARSON_P), 4))[,-1] 

#These calculations may take a few minutes to complete
OACInputData_CORR_SPEARMAN = rcorr(OACInputData_CORR, type="spearman")
OACInputData_CORR_SPEARMAN_R <- rcorr(OACInputData_CORR, type="spearman")$r 
OACInputData_CORR_SPEARMAN_P <- rcorr(OACInputData_CORR, type="spearman")$P 
OACInputData_CORR_SPEARMAN_R <- format(round(cbind(rep(-1.11, ncol(OACInputData_CORR)), OACInputData_CORR_SPEARMAN_R), 4))[,-1] 
OACInputData_CORR_SPEARMAN_P <- format(round(cbind(rep(-1.11, ncol(OACInputData_CORR)), OACInputData_CORR_SPEARMAN_P), 4))[,-1] 

##########SAVE WORKSPACE##########
gc()
save.image(paste("2011 OAC - Selecting Variables - ",ColNumber, " - Correlation Data.RData", sep=""))
##################################

#Plot Type 1
OACInputData_Correlation_Matrix_Plot <-OACInputData_P_CORR
OACInputData_Correlation_Matrix_Plot_Var_Number<-ncol(OACInputData_Correlation_Matrix_Plot)
colnames(OACInputData_Correlation_Matrix_Plot) <- c(1:OACInputData_Correlation_Matrix_Plot_Var_Number)
OACInputData_Correlation_Matrix_Plot_Cor <- cor(OACInputData_Correlation_Matrix_Plot)
OACInputData_Correlation_Matrix_Plot_Melt <- melt(OACInputData_Correlation_Matrix_Plot_Cor)
colnames(OACInputData_Correlation_Matrix_Plot_Melt) <- c("X1", "X2", "Value")
OACInputData_Correlation_Matrix_Plot_Melt$Value<-cut(OACInputData_Correlation_Matrix_Plot_Melt$Value,breaks=c(-1,-0.60,-0.20,0.20,0.6,1),include.lowest=TRUE,label=c("-1.00 to -0.60 ","-0.60 to -0.20 ","-0.20 to 0.20 ","0.20 to 0.60 ","0.60 to 1.00 "))
CorrMatrix.NoPanel <- list(theme(panel.background=element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank())) 
CorrMatrixPlot <- ggplot(OACInputData_Correlation_Matrix_Plot_Melt,aes(x = X1, y = X2))
CorrMatrixPlot <- CorrMatrixPlot + geom_tile(aes(fill=Value), colour="grey15")
CorrMatrixPlot <- CorrMatrixPlot + scale_fill_manual(values = c('-1.00 to -0.60 ' = "#B2182B", '-0.60 to -0.20 ' = "#F09B7A", '-0.20 to 0.20 ' = "#FFFFFF", '0.20 to 0.60 ' = "#87BEDA", '0.60 to 1.00 ' = "#2166AC"), name="Correlation ", drop=FALSE)
CorrMatrixPlot <- CorrMatrixPlot + scale_x_discrete(name="Variable Number", limits=c(1:OACInputData_Correlation_Matrix_Plot_Var_Number))
CorrMatrixPlot <- CorrMatrixPlot + scale_y_discrete(name="Variable Number", limits=c(1:OACInputData_Correlation_Matrix_Plot_Var_Number)) 
CorrMatrixPlot <- CorrMatrixPlot + theme(axis.title.x = element_text(size=26,colour = "black", face = "bold"),axis.text.x=element_text(angle=0, size=20))
CorrMatrixPlot <- CorrMatrixPlot + theme(axis.title.y = element_text(size=26,colour = "black", face = "bold"),axis.text.y=element_text(size=20))+ CorrMatrix.NoPanel
CorrMatrixPlot <- CorrMatrixPlot + ggtitle("Correlation Matrix ") + theme(plot.title = element_text(size = 40, colour = "black", face = "bold"))
CorrMatrixPlot <- CorrMatrixPlot + coord_equal() + theme(legend.position = "bottom") 
CorrMatrixPlot <- CorrMatrixPlot + theme(legend.title = element_text(size=26)) + theme(legend.text = element_text(size=20)) 
CorrMatrixPlot <- CorrMatrixPlot + theme(axis.title.x=element_text(hjust=0.5, vjust=0.3))
CorrMatrixPlot <- CorrMatrixPlot + theme(axis.title.y=element_text(hjust=0.5, vjust=0.4))
CorrMatrixPlot <- CorrMatrixPlot + theme(plot.title=element_text(hjust=0.5, vjust=0.9))
CorrMatrixPlot_PT1 <- CorrMatrixPlot

OACInputData_Correlation_Matrix_Plot_Cor_Row_Sums <- rowSums(OACInputData_Correlation_Matrix_Plot_Cor)
OACInputData_Correlation_Matrix_Plot_Cor_Row_Rank <-data.frame(rank(OACInputData_Correlation_Matrix_Plot_Cor_Row_Sums))
colnames(OACInputData_Correlation_Matrix_Plot_Cor_Row_Rank)<-"Rank"
OACInputData_Correlation_Matrix_Plot_Cor_Rank <-cbind(OACInputData_Correlation_Matrix_Plot_Cor, OACInputData_Correlation_Matrix_Plot_Cor_Row_Rank)
OACInputData_Correlation_Matrix_Plot_Cor_Rank_NCol <- ncol(OACInputData_Correlation_Matrix_Plot_Cor_Rank)
OACInputData_Correlation_Matrix_Plot_Cor_Order <- OACInputData_Correlation_Matrix_Plot_Cor_Rank[order(OACInputData_Correlation_Matrix_Plot_Cor_Rank[,OACInputData_Correlation_Matrix_Plot_Cor_Rank_NCol]),]
OACInputData_Correlation_Matrix_Plot_Cor_Order$Rank <-NULL
OACInputData_Correlation_Matrix_Plot_Cor_Order <- OACInputData_Correlation_Matrix_Plot_Cor_Order[,order(OACInputData_Correlation_Matrix_Plot_Cor_Order[nrow(OACInputData_Correlation_Matrix_Plot_Cor_Order),])]
OACInputData_Correlation_Matrix_Plot_Cor_Order<-as.matrix(OACInputData_Correlation_Matrix_Plot_Cor_Order)
OACInputData_Correlation_Matrix_Plot_Melt_Order <- melt(OACInputData_Correlation_Matrix_Plot_Cor_Order)
colnames(OACInputData_Correlation_Matrix_Plot_Melt_Order) <- c("X1", "X2", "Value")
OACInputData_Correlation_Matrix_Plot_Melt_Order$X1 <- factor(OACInputData_Correlation_Matrix_Plot_Melt_Order$X1, levels=unique(OACInputData_Correlation_Matrix_Plot_Melt_Order$X1),ordered=TRUE)
OACInputData_Correlation_Matrix_Plot_Melt_Order$X2 <- factor(OACInputData_Correlation_Matrix_Plot_Melt_Order$X2, levels=unique(OACInputData_Correlation_Matrix_Plot_Melt_Order$X1),ordered=TRUE)
OACInputData_Correlation_Matrix_Plot_Melt_Order$Value<-cut(OACInputData_Correlation_Matrix_Plot_Melt_Order$Value,breaks=c(-1,-0.60,-0.20,0.20,0.6,1),include.lowest=TRUE,label=c("-1.00 to -0.60 ","-0.60 to -0.20 ","-0.20 to 0.20 ","0.20 to 0.60 ","0.60 to 1.00 "))
CorrMatrix.NoPanel <- list(theme(panel.background=element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank())) 
CorrMatrixPlot_Order <- ggplot(OACInputData_Correlation_Matrix_Plot_Melt_Order,aes(x = X1, y = X2))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + geom_tile(aes(fill=Value), colour="grey15")
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + scale_fill_manual(values = c('-1.00 to -0.60 ' = "#B2182B", '-0.60 to -0.20 ' = "#F09B7A", '-0.20 to 0.20 ' = "#FFFFFF", '0.20 to 0.60 ' = "#87BEDA", '0.60 to 1.00 ' = "#2166AC"), name="Correlation ", drop=FALSE)
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + scale_x_discrete(name="Variable Number")
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + scale_y_discrete(name="Variable Number") 
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + theme(axis.title.x = element_text(size=26,colour = "black", face = "bold"),axis.text.x=element_text(angle=0, size=20))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + theme(axis.title.y = element_text(size=26,colour = "black", face = "bold"),axis.text.y=element_text(size=20))+ CorrMatrix.NoPanel
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + ggtitle("Correlation Matrix ") + theme(plot.title = element_text(size = 40, colour = "black", face = "bold"))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + coord_equal() + theme(legend.position = "bottom") 
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + theme(legend.title = element_text(size=26)) + theme(legend.text = element_text(size=20)) 
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + theme(axis.title.x=element_text(hjust=0.5, vjust=0.3))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + theme(axis.title.y=element_text(hjust=0.5, vjust=0.4))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + theme(plot.title=element_text(hjust=0.5, vjust=0.9))
CorrMatrixPlot_Order_PT1 <- CorrMatrixPlot_Order

#Plot Type 2
OACInputData_Correlation_Matrix_Plot <-OACInputData_P_CORR
OACInputData_Correlation_Matrix_Plot_Var_Number<-ncol(OACInputData_Correlation_Matrix_Plot)
colnames(OACInputData_Correlation_Matrix_Plot) <- c(1:OACInputData_Correlation_Matrix_Plot_Var_Number)
OACInputData_Correlation_Matrix_Plot_Cor <- cor(OACInputData_Correlation_Matrix_Plot)
OACInputData_Correlation_Matrix_Plot_Melt <- melt(OACInputData_Correlation_Matrix_Plot_Cor)
colnames(OACInputData_Correlation_Matrix_Plot_Melt) <- c("X1", "X2", "Value")
OACInputData_Correlation_Matrix_Plot_Melt$Value<-cut(OACInputData_Correlation_Matrix_Plot_Melt$Value,breaks=c(-1,-0.60,-0.20,0.20,0.6,1),include.lowest=TRUE,label=c("-1.00 to -0.60 ","-0.60 to -0.20 ","-0.20 to 0.20 ","0.20 to 0.60 ","0.60 to 1.00 "))
CorrMatrix.NoPanel <- list(theme(panel.background=element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank())) 
CorrMatrixPlot <- ggplot(OACInputData_Correlation_Matrix_Plot_Melt,aes(x = X1, y = X2))
CorrMatrixPlot <- CorrMatrixPlot + geom_tile(aes(fill=Value), colour="grey15")
CorrMatrixPlot <- CorrMatrixPlot + scale_fill_manual(values = c('-1.00 to -0.60 ' = "#B2182B", '-0.60 to -0.20 ' = "#F09B7A", '-0.20 to 0.20 ' = "#FFFFFF", '0.20 to 0.60 ' = "#87BEDA", '0.60 to 1.00 ' = "#2166AC"), name="Correlation ", drop=FALSE)
CorrMatrixPlot <- CorrMatrixPlot + scale_x_discrete(name="Variable Number", limits=c(1:OACInputData_Correlation_Matrix_Plot_Var_Number))
CorrMatrixPlot <- CorrMatrixPlot + scale_y_discrete(name="Variable Number", limits=c(1:OACInputData_Correlation_Matrix_Plot_Var_Number)) 
CorrMatrixPlot <- CorrMatrixPlot + theme(axis.title.x = element_text(size=26,colour = "black", face = "bold"),axis.text.x=element_text(angle=0, size=20))
CorrMatrixPlot <- CorrMatrixPlot + theme(axis.title.y = element_text(size=26,colour = "black", face = "bold"),axis.text.y=element_text(size=20))+ CorrMatrix.NoPanel
CorrMatrixPlot <- CorrMatrixPlot + ggtitle("Correlation Matrix ") + theme(plot.title = element_text(size = 40, colour = "black", face = "bold"))
CorrMatrixPlot <- CorrMatrixPlot + coord_equal() + theme(legend.position = "bottom") 
CorrMatrixPlot <- CorrMatrixPlot + theme(legend.title = element_text(size=26)) + theme(legend.text = element_text(size=20)) 
CorrMatrixPlot <- CorrMatrixPlot + theme(axis.title.x=element_text(hjust=0.5, vjust=0.3))
CorrMatrixPlot <- CorrMatrixPlot + theme(axis.title.y=element_text(hjust=0.5, vjust=0.4))
CorrMatrixPlot <- CorrMatrixPlot + theme(plot.title=element_text(hjust=0.5, vjust=0.9))
CorrMatrixPlot <- CorrMatrixPlot + theme(axis.text.x=element_text(vjust=0.5, angle=270))
CorrMatrixPlot_PT2 <- CorrMatrixPlot

OACInputData_Correlation_Matrix_Plot_Cor_Row_Sums <- rowSums(OACInputData_Correlation_Matrix_Plot_Cor)
OACInputData_Correlation_Matrix_Plot_Cor_Row_Rank <-data.frame(rank(OACInputData_Correlation_Matrix_Plot_Cor_Row_Sums))
colnames(OACInputData_Correlation_Matrix_Plot_Cor_Row_Rank)<-"Rank"
OACInputData_Correlation_Matrix_Plot_Cor_Rank <-cbind(OACInputData_Correlation_Matrix_Plot_Cor, OACInputData_Correlation_Matrix_Plot_Cor_Row_Rank)
OACInputData_Correlation_Matrix_Plot_Cor_Rank_NCol <- ncol(OACInputData_Correlation_Matrix_Plot_Cor_Rank)
OACInputData_Correlation_Matrix_Plot_Cor_Order <- OACInputData_Correlation_Matrix_Plot_Cor_Rank[order(OACInputData_Correlation_Matrix_Plot_Cor_Rank[,OACInputData_Correlation_Matrix_Plot_Cor_Rank_NCol]),]
OACInputData_Correlation_Matrix_Plot_Cor_Order$Rank <-NULL
OACInputData_Correlation_Matrix_Plot_Cor_Order <- OACInputData_Correlation_Matrix_Plot_Cor_Order[,order(OACInputData_Correlation_Matrix_Plot_Cor_Order[nrow(OACInputData_Correlation_Matrix_Plot_Cor_Order),])]
OACInputData_Correlation_Matrix_Plot_Cor_Order<-as.matrix(OACInputData_Correlation_Matrix_Plot_Cor_Order)
OACInputData_Correlation_Matrix_Plot_Melt_Order <- melt(OACInputData_Correlation_Matrix_Plot_Cor_Order)
colnames(OACInputData_Correlation_Matrix_Plot_Melt_Order) <- c("X1", "X2", "Value")
OACInputData_Correlation_Matrix_Plot_Melt_Order$X1 <- factor(OACInputData_Correlation_Matrix_Plot_Melt_Order$X1, levels=unique(OACInputData_Correlation_Matrix_Plot_Melt_Order$X1),ordered=TRUE)
OACInputData_Correlation_Matrix_Plot_Melt_Order$X2 <- factor(OACInputData_Correlation_Matrix_Plot_Melt_Order$X2, levels=unique(OACInputData_Correlation_Matrix_Plot_Melt_Order$X1),ordered=TRUE)
OACInputData_Correlation_Matrix_Plot_Melt_Order$Value<-cut(OACInputData_Correlation_Matrix_Plot_Melt_Order$Value,breaks=c(-1,-0.60,-0.20,0.20,0.6,1),include.lowest=TRUE,label=c("-1.00 to -0.60 ","-0.60 to -0.20 ","-0.20 to 0.20 ","0.20 to 0.60 ","0.60 to 1.00 "))
CorrMatrix.NoPanel <- list(theme(panel.background=element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank())) 
CorrMatrixPlot_Order <- ggplot(OACInputData_Correlation_Matrix_Plot_Melt_Order,aes(x = X1, y = X2))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + geom_tile(aes(fill=Value), colour="grey15")
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + scale_fill_manual(values = c('-1.00 to -0.60 ' = "#B2182B", '-0.60 to -0.20 ' = "#F09B7A", '-0.20 to 0.20 ' = "#FFFFFF", '0.20 to 0.60 ' = "#87BEDA", '0.60 to 1.00 ' = "#2166AC"), name="Correlation ", drop=FALSE)
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + scale_x_discrete(name="Variable Number")
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + scale_y_discrete(name="Variable Number") 
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + theme(axis.title.x = element_text(size=26,colour = "black", face = "bold"),axis.text.x=element_text(angle=0, size=20))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + theme(axis.title.y = element_text(size=26,colour = "black", face = "bold"),axis.text.y=element_text(size=20))+ CorrMatrix.NoPanel
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + ggtitle("Correlation Matrix ") + theme(plot.title = element_text(size = 40, colour = "black", face = "bold"))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + coord_equal() + theme(legend.position = "bottom") 
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + theme(legend.title = element_text(size=26)) + theme(legend.text = element_text(size=20)) 
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + theme(axis.title.x=element_text(hjust=0.5, vjust=0.3))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + theme(axis.title.y=element_text(hjust=0.5, vjust=0.4))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + theme(plot.title=element_text(hjust=0.5, vjust=0.9))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + theme(axis.text.x=element_text(vjust=0.5, angle=270))
CorrMatrixPlot_Order_PT2 <- CorrMatrixPlot_Order

#Plot Type 3
OACInputData_Correlation_Matrix_Plot <-OACInputData_P_CORR
OACInputData_Correlation_Matrix_Plot_Var_Number<-ncol(OACInputData_Correlation_Matrix_Plot)
colnames(OACInputData_Correlation_Matrix_Plot) <- c(1:OACInputData_Correlation_Matrix_Plot_Var_Number)
OACInputData_Correlation_Matrix_Plot_Cor <- cor(OACInputData_Correlation_Matrix_Plot)
OACInputData_Correlation_Matrix_Plot_Melt <- melt(OACInputData_Correlation_Matrix_Plot_Cor)
colnames(OACInputData_Correlation_Matrix_Plot_Melt) <- c("X1", "X2", "Value")
OACInputData_Correlation_Matrix_Plot_Melt$Value<-cut(OACInputData_Correlation_Matrix_Plot_Melt$Value,breaks=c(-1,-0.60,-0.20,0.20,0.6,1),include.lowest=TRUE,label=c("-1.00 to -0.60 ","-0.60 to -0.20 ","-0.20 to 0.20 ","0.20 to 0.60 ","0.60 to 1.00 "))
CorrMatrix.NoPanel <- list(theme(panel.background=element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank())) 
CorrMatrixPlot <- ggplot(OACInputData_Correlation_Matrix_Plot_Melt,aes(x = X1, y = X2))
CorrMatrixPlot <- CorrMatrixPlot + geom_tile(aes(fill=Value), colour="grey15")
CorrMatrixPlot <- CorrMatrixPlot + scale_fill_manual(values = c('-1.00 to -0.60 ' = "#B2182B", '-0.60 to -0.20 ' = "#F09B7A", '-0.20 to 0.20 ' = "#FFFFFF", '0.20 to 0.60 ' = "#87BEDA", '0.60 to 1.00 ' = "#2166AC"), name=" Correlation ", drop=FALSE)
CorrMatrixPlot <- CorrMatrixPlot + scale_x_discrete(name="Variable", limits=c(1:OACInputData_Correlation_Matrix_Plot_Var_Number),breaks = seq(1, OACInputData_Correlation_Matrix_Plot_Var_Number, 1), labels = c(1, rep("",OACInputData_Correlation_Matrix_Plot_Var_Number-2), OACInputData_Correlation_Matrix_Plot_Var_Number))
CorrMatrixPlot <- CorrMatrixPlot + scale_y_discrete(name="Variable", limits=c(1:OACInputData_Correlation_Matrix_Plot_Var_Number),breaks = seq(1, OACInputData_Correlation_Matrix_Plot_Var_Number, 1), labels = c(1, rep("",OACInputData_Correlation_Matrix_Plot_Var_Number-2), OACInputData_Correlation_Matrix_Plot_Var_Number))
CorrMatrixPlot <- CorrMatrixPlot + theme(axis.title.x = element_text(size=26,colour = "black", face = "bold"),axis.text.x=element_text(angle=0, size=20))
CorrMatrixPlot <- CorrMatrixPlot + theme(axis.title.y = element_text(size=26,colour = "black", face = "bold"),axis.text.y=element_text(size=20))+ CorrMatrix.NoPanel
CorrMatrixPlot <- CorrMatrixPlot + ggtitle("Correlation Matrix ") + theme(plot.title = element_text(size = 40, colour = "black", face = "bold"))
CorrMatrixPlot <- CorrMatrixPlot + coord_equal() + theme(legend.position = "bottom") 
CorrMatrixPlot <- CorrMatrixPlot + theme(legend.title = element_text(size=26)) + theme(legend.text = element_text(size=20)) 
CorrMatrixPlot <- CorrMatrixPlot + theme(axis.title.x=element_text(hjust=0.5, vjust=0.3))
CorrMatrixPlot <- CorrMatrixPlot + theme(axis.title.y=element_text(hjust=0.5, vjust=0.8))
CorrMatrixPlot <- CorrMatrixPlot + theme(plot.title=element_text(hjust=0.5, vjust=0.9))
CorrMatrixPlot_PT3 <- CorrMatrixPlot

OACInputData_Correlation_Matrix_Plot_Cor_Row_Sums <- rowSums(OACInputData_Correlation_Matrix_Plot_Cor)
OACInputData_Correlation_Matrix_Plot_Cor_Row_Rank <-data.frame(rank(OACInputData_Correlation_Matrix_Plot_Cor_Row_Sums))
colnames(OACInputData_Correlation_Matrix_Plot_Cor_Row_Rank)<-"Rank"
OACInputData_Correlation_Matrix_Plot_Cor_Rank <-cbind(OACInputData_Correlation_Matrix_Plot_Cor, OACInputData_Correlation_Matrix_Plot_Cor_Row_Rank)
OACInputData_Correlation_Matrix_Plot_Cor_Rank_NCol <- ncol(OACInputData_Correlation_Matrix_Plot_Cor_Rank)
OACInputData_Correlation_Matrix_Plot_Cor_Order <- OACInputData_Correlation_Matrix_Plot_Cor_Rank[order(OACInputData_Correlation_Matrix_Plot_Cor_Rank[,OACInputData_Correlation_Matrix_Plot_Cor_Rank_NCol]),]
OACInputData_Correlation_Matrix_Plot_Cor_Order$Rank <-NULL
OACInputData_Correlation_Matrix_Plot_Cor_Order <- OACInputData_Correlation_Matrix_Plot_Cor_Order[,order(OACInputData_Correlation_Matrix_Plot_Cor_Order[nrow(OACInputData_Correlation_Matrix_Plot_Cor_Order),])]
OACInputData_Correlation_Matrix_Plot_Cor_Order<-as.matrix(OACInputData_Correlation_Matrix_Plot_Cor_Order)
OACInputData_Correlation_Matrix_Plot_Melt_Order <- melt(OACInputData_Correlation_Matrix_Plot_Cor_Order)
colnames(OACInputData_Correlation_Matrix_Plot_Melt_Order) <- c("X1", "X2", "Value")
OACInputData_Correlation_Matrix_Plot_Melt_Order$X1 <- factor(OACInputData_Correlation_Matrix_Plot_Melt_Order$X1, levels=unique(OACInputData_Correlation_Matrix_Plot_Melt_Order$X1),ordered=TRUE)
OACInputData_Correlation_Matrix_Plot_Melt_Order$X2 <- factor(OACInputData_Correlation_Matrix_Plot_Melt_Order$X2, levels=unique(OACInputData_Correlation_Matrix_Plot_Melt_Order$X1),ordered=TRUE)
OACInputData_Correlation_Matrix_Plot_Melt_Order$Value<-cut(OACInputData_Correlation_Matrix_Plot_Melt_Order$Value,breaks=c(-1,-0.60,-0.20,0.20,0.6,1),include.lowest=TRUE,label=c("-1.00 to -0.60 ","-0.60 to -0.20 ","-0.20 to 0.20 ","0.20 to 0.60 ","0.60 to 1.00 "))
OACInputData_Cor_Order_Variables<-data.frame(unique(OACInputData_Correlation_Matrix_Plot_Melt_Order$X1))
colnames(OACInputData_Cor_Order_Variables) <-"Ordered Variables"
OACInputData_Cor_Order_Variables_F<-OACInputData_Cor_Order_Variables[1,]
OACInputData_Cor_Order_Variables_L<-OACInputData_Cor_Order_Variables[ColNumber,]
OACInputData_Cor_Order_Variables_F <- as.matrix(OACInputData_Cor_Order_Variables_F)
OACInputData_Cor_Order_Variables_L <- as.matrix(OACInputData_Cor_Order_Variables_L)
OACInputData_Cor_Order_Variables_First <- as.numeric(OACInputData_Cor_Order_Variables_F)
OACInputData_Cor_Order_Variables_Last <- as.numeric(OACInputData_Cor_Order_Variables_L)
CorrMatrix.NoPanel <- list(theme(panel.background=element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank())) 
CorrMatrixPlot_Order <- ggplot(OACInputData_Correlation_Matrix_Plot_Melt_Order,aes(x = X1, y = X2))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + geom_tile(aes(fill=Value), colour="grey15")
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + scale_fill_manual(values = c('-1.00 to -0.60 ' = "#B2182B", '-0.60 to -0.20 ' = "#F09B7A", '-0.20 to 0.20 ' = "#FFFFFF", '0.20 to 0.60 ' = "#87BEDA", '0.60 to 1.00 ' = "#2166AC"), name=" Correlation ", drop=FALSE)
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + scale_x_discrete(name="Variable", labels = c(OACInputData_Cor_Order_Variables_First, rep("",OACInputData_Correlation_Matrix_Plot_Var_Number-2), OACInputData_Cor_Order_Variables_Last))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + scale_y_discrete(name="Variable",labels = c(OACInputData_Cor_Order_Variables_First, rep("",OACInputData_Correlation_Matrix_Plot_Var_Number-2), OACInputData_Cor_Order_Variables_Last))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + theme(axis.title.x = element_text(size=26,colour = "black", face = "bold"),axis.text.x=element_text(angle=0, size=20))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + theme(axis.title.y = element_text(size=26,colour = "black", face = "bold"),axis.text.y=element_text(size=20))+ CorrMatrix.NoPanel
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + ggtitle("Correlation Matrix ") + theme(plot.title = element_text(size = 40, colour = "black", face = "bold"))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + coord_equal() + theme(legend.position = "bottom") 
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + theme(legend.title = element_text(size=26)) + theme(legend.text = element_text(size=20)) 
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + theme(axis.title.x=element_text(hjust=0.5, vjust=0.3))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + theme(axis.title.y=element_text(hjust=0.5, vjust=0.8))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + theme(plot.title=element_text(hjust=0.5, vjust=0.9))
CorrMatrixPlot_Order_PT3 <- CorrMatrixPlot_Order

#Plot Type 4
OACInputData_Correlation_Matrix_Plot <-OACInputData_P_CORR
OACInputData_Correlation_Matrix_Plot_Var_Number<-ncol(OACInputData_Correlation_Matrix_Plot)
colnames(OACInputData_Correlation_Matrix_Plot) <- c(1:OACInputData_Correlation_Matrix_Plot_Var_Number)
OACInputData_Correlation_Matrix_Plot_Cor <- cor(OACInputData_Correlation_Matrix_Plot)
OACInputData_Correlation_Matrix_Plot_Melt <- melt(OACInputData_Correlation_Matrix_Plot_Cor)
colnames(OACInputData_Correlation_Matrix_Plot_Melt) <- c("X1", "X2", "Value")
OACInputData_Correlation_Matrix_Plot_Melt$Value<-cut(OACInputData_Correlation_Matrix_Plot_Melt$Value,breaks=c(-1,-0.60,-0.20,0.20,0.6,1),include.lowest=TRUE,label=c("-1.00 to -0.60 ","-0.60 to -0.20 ","-0.20 to 0.20 ","0.20 to 0.60 ","0.60 to 1.00 "))
CorrMatrix.NoPanel <- list(theme(panel.background=element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank())) 
CorrMatrixPlot <- ggplot(OACInputData_Correlation_Matrix_Plot_Melt,aes(x = X1, y = X2))
CorrMatrixPlot <- CorrMatrixPlot + geom_tile(aes(fill=Value))
CorrMatrixPlot <- CorrMatrixPlot + scale_fill_manual(values = c('-1.00 to -0.60 ' = "#B2182B", '-0.60 to -0.20 ' = "#F09B7A", '-0.20 to 0.20 ' = "#FFFFFF", '0.20 to 0.60 ' = "#87BEDA", '0.60 to 1.00 ' = "#2166AC"), name=" Correlation ", drop=FALSE)
CorrMatrixPlot <- CorrMatrixPlot + scale_x_discrete(name="Variable", limits=c(1:OACInputData_Correlation_Matrix_Plot_Var_Number),breaks = seq(1, OACInputData_Correlation_Matrix_Plot_Var_Number, 1), labels = c(1, rep("",OACInputData_Correlation_Matrix_Plot_Var_Number-2), OACInputData_Correlation_Matrix_Plot_Var_Number))
CorrMatrixPlot <- CorrMatrixPlot + scale_y_discrete(name="Variable", limits=c(1:OACInputData_Correlation_Matrix_Plot_Var_Number),breaks = seq(1, OACInputData_Correlation_Matrix_Plot_Var_Number, 1), labels = c(1, rep("",OACInputData_Correlation_Matrix_Plot_Var_Number-2), OACInputData_Correlation_Matrix_Plot_Var_Number))
CorrMatrixPlot <- CorrMatrixPlot + theme(axis.title.x = element_text(size=26,colour = "black", face = "bold"),axis.text.x=element_text(angle=0, size=20))
CorrMatrixPlot <- CorrMatrixPlot + theme(axis.title.y = element_text(size=26,colour = "black", face = "bold"),axis.text.y=element_text(size=20))+ CorrMatrix.NoPanel
CorrMatrixPlot <- CorrMatrixPlot + ggtitle("Correlation Matrix ") + theme(plot.title = element_text(size = 40, colour = "black", face = "bold"))
CorrMatrixPlot <- CorrMatrixPlot + coord_equal() + theme(legend.position = "bottom") 
CorrMatrixPlot <- CorrMatrixPlot + theme(legend.title = element_text(size=26)) + theme(legend.text = element_text(size=20)) 
CorrMatrixPlot <- CorrMatrixPlot + theme(axis.title.x=element_text(hjust=0.5, vjust=0.3))
CorrMatrixPlot <- CorrMatrixPlot + theme(axis.title.y=element_text(hjust=0.5, vjust=0.8))
CorrMatrixPlot <- CorrMatrixPlot + theme(plot.title=element_text(hjust=0.5, vjust=0.9))
CorrMatrixPlot <- CorrMatrixPlot + theme(legend.background = element_rect(fill="gray85"))
CorrMatrixPlot_PT4 <- CorrMatrixPlot

OACInputData_Correlation_Matrix_Plot_Cor_Row_Sums <- rowSums(OACInputData_Correlation_Matrix_Plot_Cor)
OACInputData_Correlation_Matrix_Plot_Cor_Row_Rank <-data.frame(rank(OACInputData_Correlation_Matrix_Plot_Cor_Row_Sums))
colnames(OACInputData_Correlation_Matrix_Plot_Cor_Row_Rank)<-"Rank"
OACInputData_Correlation_Matrix_Plot_Cor_Rank <-cbind(OACInputData_Correlation_Matrix_Plot_Cor, OACInputData_Correlation_Matrix_Plot_Cor_Row_Rank)
OACInputData_Correlation_Matrix_Plot_Cor_Rank_NCol <- ncol(OACInputData_Correlation_Matrix_Plot_Cor_Rank)
OACInputData_Correlation_Matrix_Plot_Cor_Order <- OACInputData_Correlation_Matrix_Plot_Cor_Rank[order(OACInputData_Correlation_Matrix_Plot_Cor_Rank[,OACInputData_Correlation_Matrix_Plot_Cor_Rank_NCol]),]
OACInputData_Correlation_Matrix_Plot_Cor_Order$Rank <-NULL
OACInputData_Correlation_Matrix_Plot_Cor_Order <- OACInputData_Correlation_Matrix_Plot_Cor_Order[,order(OACInputData_Correlation_Matrix_Plot_Cor_Order[nrow(OACInputData_Correlation_Matrix_Plot_Cor_Order),])]
OACInputData_Correlation_Matrix_Plot_Cor_Order<-as.matrix(OACInputData_Correlation_Matrix_Plot_Cor_Order)
OACInputData_Correlation_Matrix_Plot_Melt_Order <- melt(OACInputData_Correlation_Matrix_Plot_Cor_Order)
colnames(OACInputData_Correlation_Matrix_Plot_Melt_Order) <- c("X1", "X2", "Value")
OACInputData_Correlation_Matrix_Plot_Melt_Order$X1 <- factor(OACInputData_Correlation_Matrix_Plot_Melt_Order$X1, levels=unique(OACInputData_Correlation_Matrix_Plot_Melt_Order$X1),ordered=TRUE)
OACInputData_Correlation_Matrix_Plot_Melt_Order$X2 <- factor(OACInputData_Correlation_Matrix_Plot_Melt_Order$X2, levels=unique(OACInputData_Correlation_Matrix_Plot_Melt_Order$X1),ordered=TRUE)
OACInputData_Correlation_Matrix_Plot_Melt_Order$Value<-cut(OACInputData_Correlation_Matrix_Plot_Melt_Order$Value,breaks=c(-1,-0.60,-0.20,0.20,0.6,1),include.lowest=TRUE,label=c("-1.00 to -0.60 ","-0.60 to -0.20 ","-0.20 to 0.20 ","0.20 to 0.60 ","0.60 to 1.00 "))
OACInputData_Cor_Order_Variables<-data.frame(unique(OACInputData_Correlation_Matrix_Plot_Melt_Order$X1))
colnames(OACInputData_Cor_Order_Variables) <-"Ordered Variables"
OACInputData_Cor_Order_Variables_F<-OACInputData_Cor_Order_Variables[1,]
OACInputData_Cor_Order_Variables_L<-OACInputData_Cor_Order_Variables[ColNumber,]
OACInputData_Cor_Order_Variables_F <- as.matrix(OACInputData_Cor_Order_Variables_F)
OACInputData_Cor_Order_Variables_L <- as.matrix(OACInputData_Cor_Order_Variables_L)
OACInputData_Cor_Order_Variables_First <- as.numeric(OACInputData_Cor_Order_Variables_F)
OACInputData_Cor_Order_Variables_Last <- as.numeric(OACInputData_Cor_Order_Variables_L)
CorrMatrix.NoPanel <- list(theme(panel.background=element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank())) 
CorrMatrixPlot_Order <- ggplot(OACInputData_Correlation_Matrix_Plot_Melt_Order,aes(x = X1, y = X2))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + geom_tile(aes(fill=Value))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + scale_fill_manual(values = c('-1.00 to -0.60 ' = "#B2182B", '-0.60 to -0.20 ' = "#F09B7A", '-0.20 to 0.20 ' = "#FFFFFF", '0.20 to 0.60 ' = "#87BEDA", '0.60 to 1.00 ' = "#2166AC"), name=" Correlation ", drop=FALSE)
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + scale_x_discrete(name="Variable", labels = c(OACInputData_Cor_Order_Variables_First, rep("",OACInputData_Correlation_Matrix_Plot_Var_Number-2), OACInputData_Cor_Order_Variables_Last))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + scale_y_discrete(name="Variable",labels = c(OACInputData_Cor_Order_Variables_First, rep("",OACInputData_Correlation_Matrix_Plot_Var_Number-2), OACInputData_Cor_Order_Variables_Last))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + theme(axis.title.x = element_text(size=26,colour = "black", face = "bold"),axis.text.x=element_text(angle=0, size=20))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + theme(axis.title.y = element_text(size=26,colour = "black", face = "bold"),axis.text.y=element_text(size=20))+ CorrMatrix.NoPanel
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + ggtitle("Correlation Matrix ") + theme(plot.title = element_text(size = 40, colour = "black", face = "bold"))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + coord_equal() + theme(legend.position = "bottom") 
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + theme(legend.title = element_text(size=26)) + theme(legend.text = element_text(size=20)) 
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + theme(axis.title.x=element_text(hjust=0.5, vjust=0.3))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + theme(axis.title.y=element_text(hjust=0.5, vjust=0.8))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + theme(plot.title=element_text(hjust=0.5, vjust=0.9))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + theme(legend.background = element_rect(fill="gray85"))
CorrMatrixPlot_Order_PT4 <- CorrMatrixPlot_Order

CorrMatrixPlot_OrderLookup <- merge(OACInputData_Cor_Order_Variables, OACInputData_Input_Variables, by.x = "Ordered Variables", by.y = "Variable Number", sort=F)
colnames(CorrMatrixPlot_OrderLookup) <-c("Variable Number", "Variable Name")

OACInputData_Correlation_Matrix <-OACInputData_P_CORR
OACInputData_Correlation_Matrix_Var_Number<-ncol(OACInputData_Correlation_Matrix)
colnames(OACInputData_Correlation_Matrix) <- c(1:OACInputData_Correlation_Matrix_Var_Number)
OACInputData_Correlation_Matrix_Cor <- cor(OACInputData_Correlation_Matrix)
OACInputData_Correlation_Matrix_Melt <- melt(OACInputData_Correlation_Matrix_Cor)
OACInputData_Correlation_Matrix_Melt_ID <- data.frame(1:nrow(OACInputData_Correlation_Matrix_Melt))
colnames(OACInputData_Correlation_Matrix_Melt_ID) <-"ID"
OACInputData_Correlation_Matrix_Melt_Bind <-cbind(OACInputData_Correlation_Matrix_Melt_ID, OACInputData_Correlation_Matrix_Melt)
OACInputData_Correlation_Matrix_Melt_Abs <- data.frame(abs(OACInputData_Correlation_Matrix_Melt_Bind))
OACInputData_Correlation_Matrix_Melt_Sig_Corr <- OACInputData_Correlation_Matrix_Melt_Abs[OACInputData_Correlation_Matrix_Melt_Abs$value <= 0.999999999999999,]
OACInputData_Correlation_Matrix_Melt_Sig_Corr <- OACInputData_Correlation_Matrix_Melt_Sig_Corr[OACInputData_Correlation_Matrix_Melt_Sig_Corr$value >=0.6,]
OACInputData_Correlation_Matrix_Melt_Sig_Corr <- OACInputData_Correlation_Matrix_Melt_Sig_Corr[with(OACInputData_Correlation_Matrix_Melt_Sig_Corr, order(-value)), ]
OACInputData_Correlation_Matrix_Melt_Sig_Corr$value <- NULL
colnames(OACInputData_Correlation_Matrix_Melt_Sig_Corr) <- c("ID", "X1", "X2")
colnames(OACInputData_Correlation_Matrix_Melt_Bind)<- c("ID", "X1", "X2", "value")
OACInputData_Correlation_Matrix_Melt_Sig_Corr <- merge(OACInputData_Correlation_Matrix_Melt_Sig_Corr, OACInputData_Correlation_Matrix_Melt_Bind, by= c("ID", "X1", "X2"), ,all = FALSE, sort=F)
OACInputData_Sig_Corr<- nrow(OACInputData_Correlation_Matrix_Melt_Sig_Corr)
OACInputData_Sig_Corr_Range<-data.frame(1:OACInputData_Sig_Corr)
colnames(OACInputData_Sig_Corr_Range) <-("Order")
OACInputData_Correlation_Matrix_Melt_Sig_Corr <- cbind(OACInputData_Correlation_Matrix_Melt_Sig_Corr, OACInputData_Sig_Corr_Range)
OACInputData_Correlation_Matrix_Melt_Sig_Corr_Plot<-OACInputData_Correlation_Matrix_Melt_Sig_Corr

Variable_Names<-data.frame(colnames(OACInputData_P_CORR))
Variable_Numbers<-data.frame(c(1:OACInputData_Correlation_Matrix_Var_Number))
Variable_LookupX1 <- cbind(Variable_Numbers, Variable_Names)
colnames(Variable_LookupX1) <- c("X1", "Variable 1")
Variable_LookupX2 <- cbind(Variable_Numbers, Variable_Names)
colnames(Variable_LookupX2) <- c("X2", "Variable 2")
Variable_Lookup <- merge(OACInputData_Correlation_Matrix_Melt_Sig_Corr, Variable_LookupX1, by= "X1", ,all = FALSE, sort=F)
Variable_Lookup <- merge(Variable_Lookup, Variable_LookupX2, by= "X2", ,all = FALSE, sort=F)
Sig_Variable_Lookup<-Variable_Lookup[,c("Order","Variable 1", "Variable 2", "X1", "X2", "value")]
colnames(Sig_Variable_Lookup)<- c("ID", "Variable 1 Name", "Variable 2 Name", "Variable 1 Number", "Variable 2 Number", "Correlation Value")
Sig_Variable_Lookup_All<-Sig_Variable_Lookup[with(Sig_Variable_Lookup, order(ID)), ]
Sig_Variable_Lookup_Edit<-Sig_Variable_Lookup_All[Sig_Variable_Lookup_All$ID[Sig_Variable_Lookup_All$ID %% 2 ==0],]
Sig_Variable_Lookup_Edit_ID <-1:nrow(Sig_Variable_Lookup_Edit)
Sig_Variable_Lookup_Edit$ID <-NULL
Sig_Variable_Lookup_Edit<-cbind(Sig_Variable_Lookup_Edit_ID , Sig_Variable_Lookup_Edit)
colnames(Sig_Variable_Lookup_Edit)[1] <-"ID"

OACInputData_Correlation_Matrix_Melt_Sig_Corr_Plot$value<-cut(OACInputData_Correlation_Matrix_Melt_Sig_Corr$value,breaks=c(-1,-0.9,-0.8,-0.7,0, 0.7,0.8,0.9,1),include.lowest=TRUE,label=c("-1.00 to -0.90 ","-0.90 to -0.80 ","-0.80 to -0.70 ","-0.70 to -0.60 ","0.60 to 0.70 ","0.70 to 0.80 ","0.80 to 0.90 ","0.90 to 1.00   "))

OACInputData_Correlation_Matrix_Melt_Sig_Corr_Var<-nrow(OACInputData_Correlation_Matrix_Melt_Sig_Corr)
OACInputData_Correlation_Matrix_Melt_Sig_Corr_Var_X1 <-data.frame(unique(OACInputData_Correlation_Matrix_Melt_Sig_Corr_Plot$X1))
OACInputData_Correlation_Matrix_Melt_Sig_Corr_Var_X2 <-data.frame(unique(OACInputData_Correlation_Matrix_Melt_Sig_Corr_Plot$X2))
colnames(OACInputData_Correlation_Matrix_Melt_Sig_Corr_Var_X1) <- ("UX1")
colnames(OACInputData_Correlation_Matrix_Melt_Sig_Corr_Var_X2) <- ("UX2")
OACInputData_Correlation_Matrix_Melt_Sig_Corr_Var_X1 <- OACInputData_Correlation_Matrix_Melt_Sig_Corr_Var_X1[with(OACInputData_Correlation_Matrix_Melt_Sig_Corr_Var_X1, order(UX1)), ]
OACInputData_Correlation_Matrix_Melt_Sig_Corr_Var_X2 <- OACInputData_Correlation_Matrix_Melt_Sig_Corr_Var_X2[with(OACInputData_Correlation_Matrix_Melt_Sig_Corr_Var_X2, order(UX2)), ]
Sig_Corr_Var_X1_Min <- min(OACInputData_Correlation_Matrix_Melt_Sig_Corr_Var_X1)
Sig_Corr_Var_X1_Max <- max(OACInputData_Correlation_Matrix_Melt_Sig_Corr_Var_X1)
Sig_Corr_Var_X2_Min <- min(OACInputData_Correlation_Matrix_Melt_Sig_Corr_Var_X2)
Sig_Corr_Var_X2_Max <- max(OACInputData_Correlation_Matrix_Melt_Sig_Corr_Var_X2)
Sig_Corr_Var_Count<-count(OACInputData_Correlation_Matrix_Melt_Sig_Corr_Var_X1)
Sig_Corr_Var_Count<-sum(Sig_Corr_Var_Count$freq)
Sig_Corr_Var_CNumX1<- data.frame(OACInputData_Correlation_Matrix_Melt_Sig_Corr_Var_X1)
colnames(Sig_Corr_Var_CNumX1) <- ("X1")
Sig_Corr_Var_CNumX2<- data.frame(OACInputData_Correlation_Matrix_Melt_Sig_Corr_Var_X2)
colnames(Sig_Corr_Var_CNumX2) <- ("X2")
Sig_Corr_Var_PNumZ1<- data.frame(1:Sig_Corr_Var_Count)
colnames(Sig_Corr_Var_PNumZ1) <- ("Z1")
Sig_Corr_Var_PNumZ2<- data.frame(1:Sig_Corr_Var_Count)
colnames(Sig_Corr_Var_PNumZ2) <- ("Z2")
Sig_Corr_VarX1 <-cbind(Sig_Corr_Var_CNumX1,Sig_Corr_Var_PNumZ1)
Sig_Corr_VarX2 <-cbind(Sig_Corr_Var_CNumX2,Sig_Corr_Var_PNumZ2)
OACInputData_Correlation_Matrix_Melt_Sig_Corr_Plot <- merge(OACInputData_Correlation_Matrix_Melt_Sig_Corr_Plot, Sig_Corr_VarX1, by= "X1", ,all = FALSE, sort=F)
OACInputData_Correlation_Matrix_Melt_Sig_Corr_Plot <- merge(OACInputData_Correlation_Matrix_Melt_Sig_Corr_Plot, Sig_Corr_VarX2, by= "X2", ,all = FALSE, sort=F)

CorrMatrixPlot_Sig <- ggplot(OACInputData_Correlation_Matrix_Melt_Sig_Corr_Plot,aes(x = Z1, y = Z2))
CorrMatrixPlot_Sig <- CorrMatrixPlot_Sig + geom_tile(aes(fill=value))
CorrMatrixPlot_Sig <- CorrMatrixPlot_Sig + scale_fill_manual(values = c('-1.00 to -0.90 ' = "#40004B", '-0.90 to -0.80 ' = "#762A83", '-0.80 to -0.70 ' = "#9970AB", '-0.70 to -0.60 ' = "#C2A5CF", '0.60 to 0.70 ' = "#A6DBA0", '0.70 to 0.80 ' = "#5AAE61", '0.80 to 0.90 ' = "#1B7837", '0.90 to 1.00   ' = "#00441B"), name="Correlation ", drop=FALSE)
CorrMatrixPlot_Sig <- CorrMatrixPlot_Sig + scale_x_discrete(name="Variable Number", limits=c(1:Sig_Corr_Var_Count),breaks = c(1:Sig_Corr_Var_Count), labels = OACInputData_Correlation_Matrix_Melt_Sig_Corr_Var_X1)
CorrMatrixPlot_Sig <- CorrMatrixPlot_Sig + scale_y_discrete(name="Variable Number", limits=c(1:Sig_Corr_Var_Count),breaks = c(1:Sig_Corr_Var_Count), labels = OACInputData_Correlation_Matrix_Melt_Sig_Corr_Var_X1)
CorrMatrixPlot_Sig1 <- CorrMatrixPlot_Sig + theme(axis.title.x = element_text(size=26,colour = "black", face = "bold"),axis.text.x=element_text(angle=0, size=20))
CorrMatrixPlot_Sig1 <- CorrMatrixPlot_Sig1 + theme(axis.title.y = element_text(size=26,colour = "black", face = "bold"),axis.text.y=element_text(size=20))
CorrMatrixPlot_Sig1 <- CorrMatrixPlot_Sig1 + ggtitle("Significant Correlation Matrix") + theme(plot.title = element_text(size = 40, colour = "black", face = "bold"))
CorrMatrixPlot_Sig1 <- CorrMatrixPlot_Sig1 + coord_equal() + theme(legend.position = "bottom") 
CorrMatrixPlot_Sig1 <- CorrMatrixPlot_Sig1 + theme(legend.title = element_text(size=26)) + theme(legend.text = element_text(size=15.3)) 
CorrMatrixPlot_Sig1 <- CorrMatrixPlot_Sig1 + theme(axis.title.x=element_text(hjust=0.5, vjust=0.3))
CorrMatrixPlot_Sig1 <- CorrMatrixPlot_Sig1 + theme(axis.title.y=element_text(hjust=0.5, vjust=0.4))
CorrMatrixPlot_Sig1 <- CorrMatrixPlot_Sig1 + theme(plot.title=element_text(hjust=0.5, vjust=1.6))
CorrMatrixPlot_Sig_PT1 <- CorrMatrixPlot_Sig1

CorrMatrixPlot_Sig2 <- CorrMatrixPlot_Sig + theme(axis.title.x = element_text(size=26,colour = "black", face = "bold"),axis.text.x=element_text(angle=0, size=20))
CorrMatrixPlot_Sig2 <- CorrMatrixPlot_Sig2 + theme(axis.title.y = element_text(size=26,colour = "black", face = "bold"),axis.text.y=element_text(size=20))
CorrMatrixPlot_Sig2 <- CorrMatrixPlot_Sig2 + ggtitle("Significant Correlation Matrix") + theme(plot.title = element_text(size = 40, colour = "black", face = "bold"))
CorrMatrixPlot_Sig2 <- CorrMatrixPlot_Sig2 + coord_equal() + theme(legend.position = "bottom") 
CorrMatrixPlot_Sig2 <- CorrMatrixPlot_Sig2 + theme(legend.title = element_text(size=26)) + theme(legend.text = element_text(size=15.3)) 
CorrMatrixPlot_Sig2 <- CorrMatrixPlot_Sig2 + theme(axis.title.x=element_text(hjust=0.5, vjust=0.3))
CorrMatrixPlot_Sig2 <- CorrMatrixPlot_Sig2 + theme(axis.title.y=element_text(hjust=0.5, vjust=0.4))
CorrMatrixPlot_Sig2 <- CorrMatrixPlot_Sig2 + theme(plot.title=element_text(hjust=0.5, vjust=1.6))
CorrMatrixPlot_Sig2 <- CorrMatrixPlot_Sig2 + theme(axis.text.x=element_text(vjust=0.5, angle=270))
CorrMatrixPlot_Sig_PT2 <- CorrMatrixPlot_Sig2

CorrMatrixPlot_Sig3 <- CorrMatrixPlot_Sig + theme(axis.title.x = element_text(size=26,colour = "black", face = "bold"),axis.text.x=element_text(angle=0, size=12))
CorrMatrixPlot_Sig3 <- CorrMatrixPlot_Sig3 + theme(axis.title.y = element_text(size=26,colour = "black", face = "bold"),axis.text.y=element_text(size=12))
CorrMatrixPlot_Sig3 <- CorrMatrixPlot_Sig3 + ggtitle("Significant Correlation Matrix") + theme(plot.title = element_text(size = 40, colour = "black", face = "bold"))
CorrMatrixPlot_Sig3 <- CorrMatrixPlot_Sig3 + coord_equal() + theme(legend.position = "bottom") 
CorrMatrixPlot_Sig3 <- CorrMatrixPlot_Sig3 + theme(axis.title.x=element_text(hjust=0.5, vjust=0.3))
CorrMatrixPlot_Sig3 <- CorrMatrixPlot_Sig3 + theme(axis.title.y=element_text(hjust=0.5, vjust=0.4))
CorrMatrixPlot_Sig3 <- CorrMatrixPlot_Sig3 + theme(plot.title=element_text(hjust=0.5, vjust=1.6))
CorrMatrixPlot_Sig3 <- CorrMatrixPlot_Sig3 + theme(axis.text.x=element_text(vjust=0.5, angle=270))
CorrMatrixPlot_Sig3 <- CorrMatrixPlot_Sig3 + theme(legend.text=element_text(size=15.3))
CorrMatrixPlot_Sig3 <- CorrMatrixPlot_Sig3 + theme(legend.title=element_text(size=26))
CorrMatrixPlot_Sig_PT3 <- CorrMatrixPlot_Sig3

write.table(OACInputData_Input_Variables, paste("Correlation/Correlation Matrix by Variables Order - ", OACInputData_Input_Variables[1,1], " to ", nrow(OACInputData_Input_Variables)," Variable Order.csv", sep = ""), sep = ",", row.names= FALSE, col.names = TRUE, qmethod = "double")

write.table(CorrMatrixPlot_OrderLookup, paste("Correlation/Correlation Matrix by Groupings Order - ", OACInputData_Cor_Order_Variables_First, " to ", OACInputData_Cor_Order_Variables_Last," Variable Order.csv", sep = ""), sep = ",", row.names= FALSE, col.names = TRUE, qmethod = "double")

CorrMatrixPlot_List<- list(CorrMatrixPlot_PT1, CorrMatrixPlot_PT2, CorrMatrixPlot_PT3, CorrMatrixPlot_PT4)
pdf(file = paste("Correlation/Correlation Matrix by Variables - ", ColNumber, " Variables - 4 Plot Types.pdf", sep=""), title = "Correlation Matrices by Variables", family='Courier', width=20, height=20)
CorrMatrixPlot_List
graphics.off()

CorrMatrixPlot_Order_List<- list(CorrMatrixPlot_Order_PT1, CorrMatrixPlot_Order_PT2, CorrMatrixPlot_Order_PT3, CorrMatrixPlot_Order_PT4)
pdf(file = paste("Correlation/Correlation Matrix by Groupings - ", ColNumber, " Variables - 4 Plot Types.pdf", sep=""), title = "Correlation Matrices by Groupings", family='Courier', width=20, height=20)
CorrMatrixPlot_Order_List
graphics.off()

CorrMatrixPlot_Sig_List<- list(CorrMatrixPlot_Sig_PT1, CorrMatrixPlot_Sig_PT2, CorrMatrixPlot_Sig_PT3)
pdf(file = paste("Correlation/Significant Correlation Matrix for ", ColNumber, " Variables - 3 Plot Types.pdf", sep=""), title = "Significant Correlation Matrices", family='Courier', width=20, height=20)
CorrMatrixPlot_Sig_List
graphics.off()

write.table(Sig_Variable_Lookup_All, paste("Correlation/Significant Correlation Matrix for ", ColNumber, " Variables Lookup Table - All.csv", sep = ""), sep = ",", row.names=FALSE, col.names = TRUE, qmethod = "double")

write.table(Sig_Variable_Lookup_Edit, paste("Correlation/Significant Correlation Matrix for ", ColNumber, " Variables Lookup Table - Edit.csv", sep = ""), sep = ",", row.names=FALSE, col.names = TRUE, qmethod = "double")

write.table(OACInputData_P_SEL_FINAL, paste("Correlation/Correlation and P Values - Pearson Correlation Coefficient - ", ColNumber, " Variables.csv", sep = ""), sep = ",", row.names=TRUE, col.names = NA, qmethod = "double")

write.table(OACInputData_S_SEL_FINAL, paste("Correlation/Correlation and P Values - Spearman Correlation Coefficient - ", ColNumber, " Variables.csv", sep = ""), sep = ",", row.names=TRUE, col.names = NA, qmethod = "double")

SigCorrRawS1<-Sig_Variable_Lookup_All[,2:3]
SigCorrRawS2<-Sig_Variable_Lookup_All[,6]
SigCorrRaw<-cbind(SigCorrRawS1,SigCorrRawS2)
colnames(SigCorrRaw) <-c("Var1", "Var2", "Value")
SigCorrRaw$Value<-sprintf("%.2f",SigCorrRaw$Value)
SigCorrVar2 <-data.frame(paste(SigCorrRaw$Var2, " (", SigCorrRaw$Value,")",sep=""))
colnames(SigCorrVar2) <- "High_Correlation"
SigCorrVar1<- data.frame(SigCorrRaw$Var1)
colnames(SigCorrVar1) <- "Variable"
SigCorrVariables_CorrOrder <- cbind(SigCorrVar1, SigCorrVar2)
SigCorrVariables_VarOrder<- SigCorrVariables_CorrOrder [order(SigCorrVariables_CorrOrder $High_Correlation, decreasing=FALSE),]

SigCorrVarUnique_CorrOrder <- unique(SigCorrVariables_CorrOrder$Variable)
SigCorrVarUniqueNum_CorrOrder <- length(SigCorrVarUnique_CorrOrder)
SigCorrVarUniqueOrder_CorrOrder<-sort(SigCorrVarUnique_CorrOrder, decreasing=FALSE)

SigCorrVarColMaxNum_CorrOrder <-0

for (sig in 1:SigCorrVarUniqueNum_CorrOrder)
{
SigCorrVarNum_CorrOrder<-SigCorrVarUniqueOrder_CorrOrder[sig]
SigCorrVarNumData_CorrOrder<-subset(SigCorrVariables_CorrOrder ,Variable==SigCorrVarNum_CorrOrder)
SigCorrVarNumLine_CorrOrder<-data.frame(t(SigCorrVarNumData_CorrOrder))
SigCorrVarRow_CorrOrder <- cbind(SigCorrVarNum_CorrOrder, data.frame(SigCorrVarNumLine_CorrOrder[2,]))
colnames(SigCorrVarRow_CorrOrder) <- "Variable"
row.names(SigCorrVarRow_CorrOrder) <-NULL
assign(paste("SigCorrVarRow_CorrOrder", sig, sep=""),SigCorrVarRow_CorrOrder)
SigCorrVarColNum_CorrOrder <- ncol(SigCorrVarRow_CorrOrder)
if (SigCorrVarColNum_CorrOrder > SigCorrVarColMaxNum_CorrOrder)	
		{
		SigCorrVarColMaxNum_CorrOrder <- SigCorrVarColNum_CorrOrder
		}
}

SigCorrVarColsNum_CorrOrder<- SigCorrVarColMaxNum_CorrOrder-1
SigCorrVarTable_CorrOrder <-NULL
for (corr in 1:SigCorrVarUniqueNum_CorrOrder)
{
SigCorrVarRowColNum_CorrOrder <- ncol(get(paste("SigCorrVarRow_CorrOrder",corr, sep="")))
SigCorrColBlanks_CorrOrder <- SigCorrVarColMaxNum_CorrOrder-SigCorrVarRowColNum_CorrOrder
BlankRows_CorrOrder <- matrix(c(rep.int("-",1)),nrow=1,ncol=SigCorrColBlanks_CorrOrder)
BlankRowsDF_CorrOrder <- data.frame(BlankRows_CorrOrder)
SigCorrVarRowBlanks_CorrOrder<-cbind(get(paste("SigCorrVarRow_CorrOrder",corr, sep="")), BlankRowsDF_CorrOrder)
colnames(SigCorrVarRowBlanks_CorrOrder)[2:SigCorrVarColMaxNum_CorrOrder] <-c(1:SigCorrVarColsNum_CorrOrder)
assign(paste("SigCorrVarRowBlanks_CorrOrder", corr, sep=""),SigCorrVarRowBlanks_CorrOrder)
SigCorrVarTable_CorrOrder <- rbind(SigCorrVarTable_CorrOrder, SigCorrVarRowBlanks_CorrOrder)
}

SigCorrVarUnique_VarOrder <- unique(SigCorrVariables_VarOrder$Variable)
SigCorrVarUniqueNum_VarOrder <- length(SigCorrVarUnique_VarOrder)
SigCorrVarUniqueOrder_VarOrder<-sort(SigCorrVarUnique_VarOrder, decreasing=FALSE)

SigCorrVarColMaxNum_VarOrder <-0

for (sig in 1:SigCorrVarUniqueNum_VarOrder)
{
SigCorrVarNum_VarOrder<-SigCorrVarUniqueOrder_VarOrder[sig]
SigCorrVarNumData_VarOrder<-subset(SigCorrVariables_VarOrder ,Variable==SigCorrVarNum_VarOrder)
SigCorrVarNumLine_VarOrder<-data.frame(t(SigCorrVarNumData_VarOrder))
SigCorrVarRow_VarOrder <- cbind(SigCorrVarNum_VarOrder, data.frame(SigCorrVarNumLine_VarOrder[2,]))
colnames(SigCorrVarRow_VarOrder) <- "Variable"
row.names(SigCorrVarRow_VarOrder) <-NULL
assign(paste("SigCorrVarRow_VarOrder", sig, sep=""),SigCorrVarRow_VarOrder)
SigCorrVarColNum_VarOrder <- ncol(SigCorrVarRow_VarOrder)
if (SigCorrVarColNum_VarOrder > SigCorrVarColMaxNum_VarOrder)	
		{
		SigCorrVarColMaxNum_VarOrder <- SigCorrVarColNum_VarOrder
		}
}

SigCorrVarColsNum_VarOrder<- SigCorrVarColMaxNum_VarOrder-1
SigCorrVarTable_VarOrder <-NULL
for (corr in 1:SigCorrVarUniqueNum_VarOrder)
{
SigCorrVarRowColNum_VarOrder <- ncol(get(paste("SigCorrVarRow_VarOrder",corr, sep="")))
SigCorrColBlanks_VarOrder <- SigCorrVarColMaxNum_VarOrder-SigCorrVarRowColNum_VarOrder
BlankRows_VarOrder <- matrix(c(rep.int("-",1)),nrow=1,ncol=SigCorrColBlanks_VarOrder)
BlankRowsDF_VarOrder <- data.frame(BlankRows_VarOrder)
SigCorrVarRowBlanks_VarOrder<-cbind(get(paste("SigCorrVarRow_VarOrder",corr, sep="")), BlankRowsDF_VarOrder)
colnames(SigCorrVarRowBlanks_VarOrder)[2:SigCorrVarColMaxNum_VarOrder] <-c(1:SigCorrVarColsNum_VarOrder)
assign(paste("SigCorrVarRowBlanks_VarOrder", corr, sep=""),SigCorrVarRowBlanks_VarOrder)
SigCorrVarTable_VarOrder <- rbind(SigCorrVarTable_VarOrder, SigCorrVarRowBlanks_VarOrder)
}

write.table(SigCorrVarTable_CorrOrder, paste("Correlation/Significant Correlation Lookup Table for ", ColNumber, " Variables - Correlation Order.csv", sep = ""), sep = ",", row.names=FALSE, col.names = TRUE, qmethod = "double")

write.table(SigCorrVarTable_VarOrder, paste("Correlation/Significant Correlation Lookup Table for ", ColNumber, " Variables - Variable Order.csv", sep = ""), sep = ",", row.names=FALSE, col.names = TRUE, qmethod = "double")

CorrelationEnd <- Sys.time()

CorrelationEnd-CorrelationStart

##########SAVE WORKSPACE##########
gc()
save.image(paste("2011 OAC - Selecting Variables - ",ColNumber, " - Correlation Data and Outputs.RData", sep=""))
##################################

##########PLOTS##########
StartPlots<- Sys.time()

OACInputData_Plots <-data.frame(OACInputData_Input, row.names=1)
All_OA <-rownames(OACInputData_Plots)
YMin <- min(OACInputData_Plots)
YMax <- max(OACInputData_Plots)

NoR <-nrow(OACInputData_Plots)
NoC <-ncol(OACInputData_Plots)
PlotPP<-FCol*FRow
PlotNoP<-round_any((NoC/FRow)/FCol, 1, f =ceiling)
VarList<-1:NoC

Split1<-1:PlotPP
if (PlotNoP >1)
{
if (PlotNoP >2)
{
for (j in 1:(PlotNoP-2))
{
assign(paste("Split",j+1, sep=""),(max(get(paste("Split",j, sep="")))+1):(max(get(paste("Split",j, sep="")))+PlotPP))
}
}
assign(paste("Split",PlotNoP, sep=""),(max(get(paste("Split",PlotNoP-1, sep="")))+1):max(NoC))
}

OACInputData_PlotsSplit <- OACInputData_Plots
OA <-rownames(OACInputData_PlotsSplit)

for (i in 1:PlotNoP)
{
Split<-get(paste("Split",i, sep=""))
assign(paste("OACInputData_Plots_",i, sep=""), OACInputData_PlotsSplit[Split])
}

OACInputData_PlotsSplit <- cbind(OA,OACInputData_PlotsSplit)

for (v in 1:PlotNoP)
{
OACInputData_Plots_Selection<-get(paste("OACInputData_Plots_",v, sep=""))
assign(paste("OACInputData_Plots_",v, sep=""), cbind(OA,OACInputData_Plots_Selection))
}

OACInputData_PlotsSplit <- melt(OACInputData_PlotsSplit, id.vars="OA")
colnames(OACInputData_PlotsSplit) <-c("OA", "Variable", "Value")

for (g in 1:PlotNoP)
{
OACMeltSplit<-get(paste("OACInputData_Plots_",g, sep=""))
OACMeltSplit<-melt(OACMeltSplit, id.vars="OA")
colnames(OACMeltSplit) <-c("OA", "Variable", "Value")
assign(paste("OACMeltSplit_",g, sep=""), OACMeltSplit)
}

PlotSplitMac = 0.01
if (Sys.info()['sysname']=="Windows")	
{
windowsFonts(Courier=windowsFont("TT Courier New"))
PlotSplitMac = 0.00000001
}

for (h in 1:PlotNoP)
{
OACMeltSplitPlot_OA<-get(paste("OACMeltSplit_",h, sep=""))
PlotSplit <- ggplot(OACMeltSplitPlot_OA, aes(x=OA, y=Value, group=Variable))
PlotSplit <- PlotSplit + geom_line(colour="grey20", size = PlotSplitMac, alpha = 0.5)
PlotSplit <- PlotSplit + scale_x_discrete(breaks=NULL) + theme_minimal(base_size = 12, base_family = "Courier") + theme(axis.ticks = element_blank(), axis.text.y = element_blank()) + theme(axis.title.y = element_blank(), axis.title.x = element_blank(),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank()) + scale_y_continuous(limits=c(YMin, YMax))
PlotSplit <- PlotSplit + facet_wrap( ~ Variable, ncol=FCol, nrow=FRow) + theme(strip.text.x = element_text(size=14)) + theme(aspect.ratio = 9/16)
assign(paste("PlotSplit_",h, sep=""), PlotSplit)
}

for (e in 1:PlotNoP)
{
OACMeltSplitPlot_Hist<-get(paste("OACMeltSplit_",e, sep=""))
PlotHist <- ggplot(OACMeltSplitPlot_Hist, aes(x=Value))
PlotHist<- PlotHist + geom_density() + theme_minimal(base_size = 12, base_family = "Courier") 
PlotHist<- PlotHist + theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(),axis.title.x = element_blank(),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
PlotHist<- PlotHist + facet_wrap(~ Variable, ncol=FCol, nrow=FRow, scales="free") + theme(strip.text.x = element_text(size=14)) + theme(aspect.ratio = 9/16)
assign(paste("PlotHist_",e, sep=""), PlotHist)
}

for (d in 1:PlotNoP)
{
PlotSplit_ggsave<-get(paste("PlotSplit_",d, sep=""))
ggsave(PlotSplit_ggsave,filename=paste("Distribution Plots/Output Area Distribution Plot ",d," of ", PlotNoP,".png",sep=""), width=30, height=20, units = c("cm"), dpi=300) 
}

for (t in 1:PlotNoP)
{
PlotHist_ggsave<-get(paste("PlotHist_",t, sep=""))
ggsave(PlotHist_ggsave,filename=paste("Distribution Plots/Density Curve Distribution Plot ",t," of ", PlotNoP,".png",sep=""), width=30, height=20, units = c("cm"), dpi=300) 
}

EndPlots <- Sys.time()


##########SAVE WORKSPACE##########
gc()
save.image(paste("2011 OAC - Selecting Variables - ",ColNumber, " - All Correlation abd Plot Data.RData", sep=""))
##################################


##########CLUSTERING WITH MISSING VARIABLES##########

MissingStart <- Sys.time()

OACInputData<-OACInputData_Input[,2:ncol(OACInputData_Input)]

OACInputData_Var_ID <- data.frame (1:nrow(OACInputData),OACInputData)
colnames(OACInputData_Var_ID) <- c("Variable_ID",names(OACInputData))

VNum <- ncol(OACInputData)

MaxCN<-max(CN)
MinCN<-min(CN)

for (r in CN)
{

OACInputData_Best_WSS_N <-1000000000000
OACInputData_Best_WSS_A <-c(1:5)
OACInputData_Best_Miss_Var_A <-c(1:r)
count=1
VarNumLoop=1

	for (u in names(OACInputData))
	{
	OACInputData_Sel<- OACInputData[, names(OACInputData)!= u]

		for (i in 1:KM)
		{
		Cluster_OACInputData <- kmeans(OACInputData_Sel, r)
		OACInputData_Best_WSS_B <- mean(Cluster_OACInputData$withinss)
		OACInputData_Best_Var_Miss_C <-data.frame (i, OACInputData_Best_WSS_B, Cluster_OACInputData$tot.withinss, Cluster_OACInputData$betweenss, Cluster_OACInputData$totss)
		OACInputData_Best_Var_Miss_B <-rbind (OACInputData_Best_WSS_A, OACInputData_Best_Var_Miss_C)
			if (OACInputData_Best_WSS_B < OACInputData_Best_WSS_N)	
			{
			BestWSS <- Cluster_OACInputData
			OACInputData_Best_WSS_N <- OACInputData_Best_WSS_B
			}
		cat("Clustering with Missing Variable", VarNumLoop, "of", VNum, "-", i, "times of", max(1:KM), "- for", r, "of", MaxCN, "Clusters \r")
		flush.console()
		}

	OACInputData_Best_Var_Miss_B <- OACInputData_Best_Var_Miss_B[2:nrow(OACInputData_Best_Var_Miss_B),]

	Best_Cluster_Group <-BestWSS$cluster

	OACInputData_Best_Var_Miss_A <- data.frame(OACInputData_Sel, Best_Cluster_Group)
	
		if (RQMISSOUT == "YES")
		{
		write.table(OACInputData_Best_Var_Miss_A, paste("Missing Variables/", count, "_OACInputData_KMeans_Missing_", u,"_",r, "_Clusters.csv", sep = ""), sep = "," , col.names = NA, qmethod = "double")
		}
		if (RQCLUSUP == "YES")
		{
		cat(paste("Current_progress:Variable_", count, "_of_", VNum,"_for_",r, "_of_", MaxCN, "_Clusters", sep = ""),file="Missing Variables/Cluster Update/Cluster Update.txt")
		}

	OACInputData_Best_Miss_Var_B<- data.frame (count, u, i, OACInputData_Best_WSS_N, Cluster_OACInputData$tot.withinss, Cluster_OACInputData$betweenss, Cluster_OACInputData$totss)
	OACInputData_Best_WSS_N <-1000000000000
	OACInputData_Best_Miss_Var_A <- rbind(OACInputData_Best_Miss_Var_A, OACInputData_Best_Miss_Var_B)
	count=count+1
	VarNumLoop=VarNumLoop+1
	}

OACInputData_Best_Miss_Var_A <- OACInputData_Best_Miss_Var_A[2:nrow(OACInputData_Best_Miss_Var_A),]

colnames(OACInputData_Best_Miss_Var_A) <- c( 'Variable_Number', 'Variable_Missing', 'KMeans_Runs', 'Lowest_Within_Cluster_Sum_of_Squares', 'Total_Within_Cluster_Sum_of_Squares', 'Total_Between_Cluster_Sum_of_Squares', 'Total_Within_&_Between_Cluster_Sum_of_Squares')

assign(paste("OACInputData_Best_Miss_Var_A_", r, sep=""), OACInputData_Best_Miss_Var_A)
}

VarNumbers <-OACInputData_Best_Miss_Var_A[,1]
VarNumbers <- data.frame(VarNumbers)
colnames(VarNumbers) <- "Missing_Variable"
LowestWCSS <- VarNumbers

for (r in CN)
{
Data<-c()
Data<-data.frame(get(paste("OACInputData_Best_Miss_Var_A_",r, sep="")))
Data_Low <- Data[,4]
Data_Low <- data.frame(Data_Low)
colnames(Data_Low) <- paste(r)
LowestWCSS <-cbind(LowestWCSS, Data_Low)
}

LowestWCSS$Mean=rowMeans(LowestWCSS[,c(2:ncol(LowestWCSS))], na.rm=TRUE)

##########SAVE WORKSPACE##########
gc()
save.image(paste("2011 OAC - Selecting Variables - ",ColNumber, " - All Correlation and Missing Variable Clustering Data.RData", sep=""))
##################################

Badger=0

for(iop in 1:ColNumber)
{
Fox<-try(get(paste('G',iop, sep="")),silent=TRUE)
	if (is.integer(Fox)=="TRUE")
	{
	Badger=Badger+1
	}
rm(Fox)
}

Pig <- Badger-1

G1_Max<-max(get(paste("G1")))
assign(paste("G1_MissingMax", sep=""),get(paste("G1_Max"))+0.49)

for(lamb in 2:Pig)
{
assign(paste("G", lamb, "_Max", sep=""),max(get(paste("G",lamb,sep=""))))
assign(paste("G", lamb, "_MissingMax", sep=""),get(paste("G", lamb, "_Max", sep=""))+0.49)
}

for(lamb in 2:Pig)
{
assign(paste("G", lamb, "_Min", sep=""),min(get(paste("G",lamb,sep=""))))
assign(paste("G", lamb, "_MissingMin", sep=""),get(paste("G", lamb, "_Min", sep=""))-0.49)
}

assign(paste("G", Badger, "_Min", sep=""),min(get(paste("G",Badger,sep=""))))
assign(paste("G", Badger, "_MissingMin", sep=""),get(paste("G", Badger, "_Min", sep=""))-0.49)

NewRowAll<-c()

for(bull in 1:Pig)
{

RowMax<-get(paste("G",bull,"_Max", sep=""))
RowMin<-get(paste("G",bull+1,"_Min", sep=""))

RowMissingMax<-get(paste("G",bull,"_MissingMax", sep=""))
RowMissingMin<-get(paste("G",bull+1,"_MissingMin", sep=""))

NewRow1 <-(LowestWCSS[RowMin,]-LowestWCSS[RowMax,])/2
NewRow2 <-LowestWCSS[RowMax,]+NewRow1

NewRow2[1]<-NULL

NewRowMax<-cbind(RowMissingMax,NewRow2)
rownames(NewRowMax)<- NewRowMax[,1]
colnames(NewRowMax)[1]<-"Missing_Variable"

NewRowMin<-cbind(RowMissingMin,NewRow2)
rownames(NewRowMin)<- NewRowMin[,1]
colnames(NewRowMin)[1]<-"Missing_Variable"

NewRowAll <-rbind(NewRowAll,NewRowMax,NewRowMin)
}

LowestWCSSAll<- rbind(LowestWCSS, NewRowAll)

LowestWCSSAll<- LowestWCSSAll[with(LowestWCSSAll, order(Missing_Variable)), ]

pdata<-melt(LowestWCSSAll, id.vars="Missing_Variable")
jdata <-pdata

LL<-list(c(unlist(G1),G1_MissingMax))
G1_New<- as.numeric(unlist(LL))

for(piglet in 2:Pig)
{
G_List <-get(paste("G",piglet, sep=""))
ListMin <-get(paste("G",piglet,"_MissingMin", sep=""))
ListMax <-get(paste("G",piglet,"_MissingMax", sep=""))

LL<-list(c(unlist(G_List),ListMin, ListMax))
assign(paste("G",piglet,"_New", sep=""),as.numeric(unlist(LL)))
}

G_List <-get(paste("G",Badger, sep=""))
ListMin <-get(paste("G",Badger,"_MissingMin", sep=""))

LL<-list(c(unlist(G_List),ListMin))
assign(paste("G",Badger,"_New", sep=""),as.numeric(unlist(LL)))

pdata_missing<-as.data.frame(pdata$Missing_Variable)

for(cow in 1:Badger)
{
groupname<-get(paste("N",cow,sep=""))
grouprange<-get(paste("G",cow,"_New",sep=""))
pdata_missing <- colwise(recode)(pdata_missing,recodes="grouprange=groupname")
}

colnames(pdata_missing)<-c("group")
pdata<-cbind(pdata,pdata_missing)

pdata_m<-subset(pdata,variable =='Mean')

VarArea <- data.frame(unique(pdata$Missing_Variable))
VarFinalArea <-c()

for (g in 1:max(VarArea))
{
VarNum<-pdata[which(pdata$Missing_Variable==g),1:c(ncol(pdata))]
VarMin <- min(VarNum$value)
VarMax <- max(VarNum$value)
VarMinMax <- cbind(g, VarMin, VarMax)
VarFinalArea <- rbind(VarFinalArea,VarMinMax)
}

for(g in 1:max(nrow(NewRowAll)))
{
frog<-NewRowAll[g,1]
VarNum<-pdata[which(pdata$Missing_Variable==frog),1:c(ncol(pdata))]
VarMin <- min(VarNum$value)
VarMax <- max(VarNum$value)
VarMinMax <- cbind(frog, VarMin, VarMax)
VarFinalArea <- rbind(VarFinalArea,VarMinMax)
}

VarFinalArea <-as.data.frame(VarFinalArea)
VarFinalArea<- VarFinalArea[with(VarFinalArea, order(g)), ]

VarFinal <- VarFinalArea
colnames(VarFinal)<-c("Variable", "Min", "Max")

varfinal_missing<-as.data.frame(VarFinal$Variable)

for(cow in 1:Badger)
{
groupname<-get(paste("N",cow,sep=""))
grouprange<-get(paste("G",cow,"_New",sep=""))
varfinal_missing <- colwise(recode)(varfinal_missing,recodes="grouprange=groupname")
}

colnames(varfinal_missing)<-c("group")
VarFinal<-cbind(VarFinal,varfinal_missing)

VarArea <- data.frame(unique(jdata$Missing_Variable))
NumV <- max(VarArea)
VarFinalMin <-c()
for (g in 1:NumV)
{
VarNum<-jdata[which(jdata$Missing_Variable==g),1:c(ncol(jdata))]
VarN <- g
VarMin <- min(VarNum$value)
VarMinN <- cbind(VarN, VarMin)
VarFinalMin <- rbind(VarFinalMin,VarMinN)
}

for(g in 1:max(nrow(NewRowAll)))
{
frog<-NewRowAll[g,1]
VarNum<-jdata[which(jdata$Missing_Variable==frog),1:c(ncol(jdata))]
VarN <- g
VarMin <- min(VarNum$value)
VarMinN <- cbind(frog, VarMin)
VarFinalMin <- rbind(VarFinalMin,VarMinN)
}

VarFinalMin <-as.data.frame(VarFinalMin)
VarFinalMin<- VarFinalMin[with(VarFinalMin, order(VarN)), ]
VFinalMin <- data.frame(VarFinalMin)
colnames(VFinalMin)<-c("Variable", "Value")
VFinalMin_Max<-max(VFinalMin$Value)

VarFinalMax <-c()
for (g in 1:NumV)
{
VarNum<-jdata[which(jdata$Missing_Variable==g),1:c(ncol(jdata))]
VarN <- g
VarMax <- max(VarNum$value)
VarMaxN <- cbind(VarN, VarMax)
VarFinalMax <- rbind(VarFinalMax,VarMaxN)
}

for(g in 1:max(nrow(NewRowAll)))
{
frog<-NewRowAll[g,1]
VarNum<-jdata[which(jdata$Missing_Variable==frog),1:c(ncol(jdata))]
VarN <- g
VarMax <- max(VarNum$value)
VarMaxN <- cbind(frog, VarMax)
VarFinalMax <- rbind(VarFinalMax,VarMaxN)
}

VarFinalMax <-as.data.frame(VarFinalMax)
VarFinalMax<- VarFinalMax[with(VarFinalMax, order(VarN)), ]
VFinalMax <- data.frame(VarFinalMax)
colnames(VFinalMax)<-c("Variable", "Value")
VFinalMax_Min<-min(VFinalMax$Value)

VFinal<-rbind(VFinalMin, VFinalMax)

vfinal_missing<-as.data.frame(VFinal$Variable)

for(cow in 1:Badger)
{
groupname<-get(paste("N",cow,sep=""))
grouprange<-get(paste("G",cow,"_New",sep=""))
vfinal_missing <- colwise(recode)(vfinal_missing,recodes="grouprange=groupname")
}

colnames(vfinal_missing)<-c("Group")
VFinal<-cbind(VFinal,vfinal_missing)

VLarge <- VFinal[which(VFinal$Value >= VFinalMax_Min-1),]
VLarge <- VLarge[order(-VLarge$Variable, VLarge$Group),]
VSmall <- VFinal[which(VFinal$Value <= VFinalMin_Max+1),]
VSmall <- VSmall[order(VSmall$Variable, VSmall$Group),]
VFinalAll <- rbind(VSmall, VLarge)

AMin <-min(jdata$value)
AMax <-max(jdata$value)
AMinRound<-round_any(AMin, 1)
AMaxRound<-round_any(AMax, 1)
AMaxAxis<-(AMax-AMin)/5
AMaxChar<-nchar(sprintf("%.0f",AMaxAxis))
AMaxCharMaxus1<-AMaxChar-1
AMaxFinalAxis <-1
for (j1 in 1:AMaxCharMaxus1)
{
AMaxFinalAxis<-paste(AMaxFinalAxis, "0", sep = "")
}
AMaxFinalAxis<-as.numeric(AMaxFinalAxis)
AMaxFinal<-round_any(AMaxAxis, AMaxFinalAxis)
AMinValue<-round_any(AMin,AMaxFinal, f =floor)
AMaxValue<-round_any(AMax,AMaxFinal, f =ceiling)

MMin <-min(pdata_m$value)
MMax <-max(pdata_m$value)
MMinRound<-round_any(MMin, 1)
MMaxRound<-round_any(MMax, 1)
MMaxAxis<-(MMax-MMin)/5
MMaxChar<-nchar(sprintf("%.0f",MMaxAxis))
MMaxCharMaxus1<-MMaxChar-1
MMaxFinalAxis <-1
for (j1 in 1:MMaxCharMaxus1)
{
MMaxFinalAxis<-paste(MMaxFinalAxis, "0", sep = "")
}
MMaxFinalAxis<-as.numeric(MMaxFinalAxis)
MMaxFinal<-round_any(MMaxAxis, MMaxFinalAxis)
MMinValue<-round_any(MMin,MMaxFinal, f =floor)
MMaxValue<-round_any(MMax,MMaxFinal, f =ceiling)

VFinal_Mean <- VFinalAll
VFinal_MeanRows <-nrow(VFinal_Mean)/2
VFinal_MeanMatchMax <- count(VFinal_Mean$Value[(VFinal_Mean$Value >=MMax) == "TRUE"])
VFinal_MeanMatchMin <- count(VFinal_Mean$Value[(VFinal_Mean$Value <=MMin) == "TRUE"])
VFinal_MatchMax <-sum(VFinal_MeanMatchMax$freq)
VFinal_MatchMin <-sum(VFinal_MeanMatchMin$freq)

if (VFinal_MatchMax == VFinal_MeanRows & VFinal_MatchMin == VFinal_MeanRows)
{
VFinal_Mean$Value[(VFinal_Mean$Value >=MMax) == "TRUE"] <- MMaxValue
VFinal_Mean$Value[(VFinal_Mean$Value <=MMin) == "TRUE"] <- MMinValue
}else
{
VFM<-VFinal_Mean
VFMTop <- VFM[seq(1, unique(VFinal_MeanRows), 1), ]
VFMBottom <- VFM[seq(unique(VFinal_MeanRows)+1, nrow(VFM), 1), ]
VFMTop$Value <- MMinValue
VFMBottom$Value <- MMaxValue
VFM <-rbind(VFMTop,VFMBottom)
VFinal_Mean <-VFM
}

VFinal_MeanAll <-VFinal_Mean

pdata<-melt(LowestWCSS, id.vars="Missing_Variable")
jdata <-pdata

pdata_missing<-as.data.frame(pdata$Missing_Variable)

for(cow in 1:Badger)
{
groupname<-get(paste("N",cow,sep=""))
grouprange<-get(paste("G",cow,sep=""))
pdata_missing <- colwise(recode)(pdata_missing,recodes="grouprange=groupname")
}

colnames(pdata_missing)<-c("group")
pdata<-cbind(pdata,pdata_missing)

pdata_m<-subset(pdata,variable =='Mean')

VarArea <- data.frame(unique(pdata$Missing_Variable))
VarFinalArea <-c()
for (g in 1:max(VarArea))
{
VarNum<-pdata[which(pdata$Missing_Variable==g),1:c(ncol(pdata))]
VarMin <- min(VarNum$value)
VarMax <- max(VarNum$value)
VarMinMax <- cbind(VarMin, VarMax)
VarFinalArea <- rbind(VarFinalArea,VarMinMax)
}
VarFinal <- data.frame(cbind(VarArea,VarFinalArea))
colnames(VarFinal)<-c("Variable", "Min", "Max")

varfinal_missing<-as.data.frame(VarFinal$Variable)

for(cow in 1:Badger)
{
groupname<-get(paste("N",cow,sep=""))
grouprange<-get(paste("G",cow,sep=""))
varfinal_missing <- colwise(recode)(varfinal_missing,recodes="grouprange=groupname")
}

colnames(varfinal_missing)<-c("group")
VarFinal<-cbind(VarFinal,varfinal_missing)


VarArea <- data.frame(unique(jdata$Missing_Variable))
NumV <- max(VarArea)
VarFinalMin <-c()
for (g in 1:NumV)
{
VarNum<-jdata[which(jdata$Missing_Variable==g),1:c(ncol(jdata))]
VarN <- g
VarMin <- min(VarNum$value)
VarMinN <- cbind(VarN, VarMin)
VarFinalMin <- rbind(VarFinalMin,VarMinN)
}
VFinalMin <- data.frame(VarFinalMin)
colnames(VFinalMin)<-c("Variable", "Value")
VFinalMin_Max<-max(VFinalMin$Value)

VarFinalMax <-c()
for (g in 1:NumV)
{
VarNum<-jdata[which(jdata$Missing_Variable==g),1:c(ncol(jdata))]
VarN <- g
VarMax <- max(VarNum$value)
VarMaxN <- cbind(VarN, VarMax)
VarFinalMax <- rbind(VarFinalMax,VarMaxN)
}
VFinalMax <- data.frame(VarFinalMax)
colnames(VFinalMax)<-c("Variable", "Value")
VFinalMax_Min<-min(VFinalMax$Value)

VFinal<-rbind(VFinalMin, VFinalMax)

vfinal_missing<-as.data.frame(VFinal$Variable)

for(cow in 1:Badger)
{
groupname<-get(paste("N",cow,sep=""))
grouprange<-get(paste("G",cow,sep=""))
vfinal_missing <- colwise(recode)(vfinal_missing,recodes="grouprange=groupname")
}

colnames(vfinal_missing)<-c("Group")
VFinal<-cbind(VFinal,vfinal_missing)

if (RQCUSTOMV == "YES")
{
VFinal <- vdata
}

VLarge <- VFinal[which(VFinal$Value >= VFinalMax_Min-1),]
VLarge <- VLarge[order(-VLarge$Variable, VLarge$Group),]
VSmall <- VFinal[which(VFinal$Value <= VFinalMin_Max+1),]
VSmall <- VSmall[order(VSmall$Variable, VSmall$Group),]
VFinal <- rbind(VSmall, VLarge)

AMin <-min(jdata$value)
AMax <-max(jdata$value)
AMinRound<-round_any(AMin, 1)
AMaxRound<-round_any(AMax, 1)
AMaxAxis<-(AMax-AMin)/5
AMaxChar<-nchar(sprintf("%.0f",AMaxAxis))
AMaxCharMaxus1<-AMaxChar-1
AMaxFinalAxis <-1
for (j1 in 1:AMaxCharMaxus1)
{
AMaxFinalAxis<-paste(AMaxFinalAxis, "0", sep = "")
}
AMaxFinalAxis<-as.numeric(AMaxFinalAxis)
AMaxFinal<-round_any(AMaxAxis, AMaxFinalAxis)
AMinValue<-round_any(AMin,AMaxFinal, f =floor)
AMaxValue<-round_any(AMax,AMaxFinal, f =ceiling)

MMin <-min(pdata_m$value)
MMax <-max(pdata_m$value)
MMinRound<-round_any(MMin, 1)
MMaxRound<-round_any(MMax, 1)
MMaxAxis<-(MMax-MMin)/5
MMaxChar<-nchar(sprintf("%.0f",MMaxAxis))
MMaxCharMaxus1<-MMaxChar-1
MMaxFinalAxis <-1
for (j1 in 1:MMaxCharMaxus1)
{
MMaxFinalAxis<-paste(MMaxFinalAxis, "0", sep = "")
}
MMaxFinalAxis<-as.numeric(MMaxFinalAxis)
MMaxFinal<-round_any(MMaxAxis, MMaxFinalAxis)
MMinValue<-round_any(MMin,MMaxFinal, f =floor)
MMaxValue<-round_any(MMax,MMaxFinal, f =ceiling)

VFinal_Mean <- VFinal
VFinal_MeanRows <-nrow(VFinal_Mean)/2
VFinal_MeanMatchMax <- count(VFinal_Mean$Value[(VFinal_Mean$Value >=MMax) == "TRUE"])
VFinal_MeanMatchMin <- count(VFinal_Mean$Value[(VFinal_Mean$Value <=MMin) == "TRUE"])
VFinal_MatchMax <-sum(VFinal_MeanMatchMax$freq)
VFinal_MatchMin <-sum(VFinal_MeanMatchMin$freq)

if (VFinal_MatchMax == VFinal_MeanRows & VFinal_MatchMin == VFinal_MeanRows)
{
VFinal_Mean$Value[(VFinal_Mean$Value >=MMax) == "TRUE"] <- MMaxValue
VFinal_Mean$Value[(VFinal_Mean$Value <=MMin) == "TRUE"] <- MMinValue
}else
{
VFM<-VFinal_Mean
VFMTop <- VFM[seq(1, unique(VFinal_MeanRows), 1), ]
VFMBottom <- VFM[seq(unique(VFinal_MeanRows)+1, nrow(VFM), 1), ]
VFMTop$Value <- MMinValue
VFMBottom$Value <- MMaxValue
VFM <-rbind(VFMTop,VFMBottom)
VFinal_Mean <-VFM
}

#Outputs for individual Variable Groups
UGroups<-unique(pdata$group)
UGroupsFreq<-count(UGroups)
UGroups<-sum(UGroupsFreq$freq)

for (o in 1:UGroups)
{
USubset<-unique(pdata$group)[o]
assign(paste("pdata_subset",o, sep=""), pdata[which(pdata$group==USubset),1:c(ncol(pdata))])
pdata_subsetmeanorder<-assign(paste("pdata_subsetmean",o, sep=""), subset(get(paste("pdata_subset",o, sep="")),variable=='Mean'))

pdata_subsethightolow<-pdata_subsetmeanorder[with(pdata_subsetmeanorder, order(-value)),]
pdata_subsetlowtohigh<-pdata_subsetmeanorder[with(pdata_subsetmeanorder, order(value)),]
pdata_subsethightolow$variable<-NULL
pdata_subsetlowtohigh$variable<-NULL
pdata_subsethightolow<-pdata_subsethightolow[c(3,1,2)]
pdata_subsetlowtohigh<-pdata_subsetlowtohigh[c(3,1,2)]
colnames(pdata_subsethightolow) <- c("Variable Domain", "Missing Variable", "WCSSValue")
colnames(pdata_subsetlowtohigh) <- c("Variable Domain", "Missing Variable", "WCSSValue")
assign(paste("pdata_subsethightolow",o, sep=""),pdata_subsethightolow)
assign(paste("pdata_subsetlowtohigh",o, sep=""),pdata_subsetlowtohigh)
}

Missing_Variables_WCSS_High_to_Low_Means_by_Group <-c()
Missing_Variables_WCSS_Low_to_High_Means_by_Group <-c()
for (t in 1:UGroups)
{
Missing_Variables_WCSS_High_to_Low_Means_by_Group<-rbind(Missing_Variables_WCSS_High_to_Low_Means_by_Group,get(paste("pdata_subsethightolow",t, sep="")))
Missing_Variables_WCSS_Low_to_High_Means_by_Group<-rbind(Missing_Variables_WCSS_Low_to_High_Means_by_Group,get(paste("pdata_subsetlowtohigh",t, sep="")))
}

Missing_Variables_WCSS_Low_to_High_Means_All_Variables<-Missing_Variables_WCSS_Low_to_High_Means_by_Group[with(Missing_Variables_WCSS_Low_to_High_Means_by_Group, order(WCSSValue)),]
Missing_Variables_WCSS_High_to_Low_Means_All_Variables<-Missing_Variables_WCSS_High_to_Low_Means_by_Group[with(Missing_Variables_WCSS_High_to_Low_Means_by_Group, order(-WCSSValue)),]
colnames(Missing_Variables_WCSS_Low_to_High_Means_by_Group) <- c("Variable Domain", "Missing Variable", "WCSS Value")
colnames(Missing_Variables_WCSS_High_to_Low_Means_by_Group) <- c("Variable Domain", "Missing Variable", "WCSS Value")
colnames(Missing_Variables_WCSS_Low_to_High_Means_All_Variables) <- c("Variable Domain", "Missing Variable", "WCSS Value")
colnames(Missing_Variables_WCSS_High_to_Low_Means_All_Variables) <- c("Variable Domain", "Missing Variable", "WCSS Value")

write.table(Missing_Variables_WCSS_Low_to_High_Means_by_Group, paste("Missing Variables/Summary Data/Missing Variables WCSS Low to High Means by Domain - ", MinCN, " to ", MaxCN," Clusters.csv", sep = ""), sep = ",", row.names= FALSE, col.names = TRUE, qmethod = "double")
write.table(Missing_Variables_WCSS_High_to_Low_Means_by_Group, paste("Missing Variables/Summary Data/Missing Variables WCSS High to Low Means by Domain - ", MinCN, " to ", MaxCN," Clusters.csv", sep = ""), sep = ",", row.names= FALSE, col.names = TRUE, qmethod = "double")
write.table(Missing_Variables_WCSS_Low_to_High_Means_All_Variables, paste("Missing Variables/Summary Data/Missing Variables WCSS Low to High Means for All Variables - ", MinCN, " to ", MaxCN," Clusters.csv", sep = ""), sep = ",", row.names= FALSE, col.names = TRUE, qmethod = "double")
write.table(Missing_Variables_WCSS_High_to_Low_Means_All_Variables, paste("Missing Variables/Summary Data/Missing Variables WCSS High to Low Means for All Variables - ", MinCN, " to ", MaxCN," Clusters.csv", sep = ""), sep = ",", row.names= FALSE, col.names = TRUE, qmethod = "double")

pdata_finalplotList <-list()

for (d in 1:UGroups)
{
toplot<-get(paste("pdata_subset",d, sep=""))
MinAxis <- min(toplot$Missing_Variable)
MaxAxis <- max(toplot$Missing_Variable)
AxisLength <- (MaxAxis-MinAxis)+1
newvalues <-c(1:AxisLength)
curvalues <- c(MinAxis:MaxAxis)
toplot$newVar <- newvalues[match(toplot$Missing_Variable, curvalues)]
toplothead <-unique(toplot$group)

cplot <-assign(paste("pdata_plot",d, sep=""), ggplot(toplot, aes(variable, value, colour = group)) +scale_x_discrete(name="Missing Variable Number", limits=c(1:AxisLength),labels=c(MinAxis:MaxAxis)) + scale_y_continuous(name="Within Cluster Sum of Squares", limits=c(AMinValue, AMaxValue), breaks = seq(AMinValue, AMaxValue, AMaxFinal)))

cplot <- cplot +geom_line(data=subset(toplot,variable =='Mean'), alpha=1, aes(newVar, value, size= ''), colour="gray20", inherit.aes = FALSE)+ guides (size = guide_legend(title = "Mean", override.aes = list(size = 3, colour = "gray20")))

cplot <- cplot +geom_line(colour="grey20", data=subset(toplot,variable =='Mean'),size = 3 , alpha=1, aes(newVar, value), inherit.aes = FALSE)

cplot <- cplot +geom_line(data=subset(toplot,variable!='Mean'),size=0.5, alpha=0.5, aes(newVar, value, shape=variable), colour='black', inherit.aes = FALSE) 

cplot <- cplot +geom_point(data=subset(toplot,variable!='Mean'), aes(newVar, value, shape = variable, col=group), size = 5.0, alpha=1.0, inherit.aes = FALSE) + labs (shape = "Clusters") + guides(scale_alpha(guide='none')) + guides(colour = "none") + guides(shape = guide_legend(override.aes = list(size = 10))) + scale_colour_manual(values = c('Demographic' = C1, 'Household Composition' = C2, 'Housing' = C3, 'Socio-Economic' = C4, 'Employment' = C5))  + theme(legend.key.height=unit(6,"line"))+ theme(legend.key.height=unit(6,"line"))

cplot <- cplot + ggtitle(paste("Clusters with Missing Variables -", toplothead, "Domain", sep=" ")) + theme(plot.title = element_text(size = 40, colour = "black", face = "bold"))+ theme(legend.title = element_text(size=26)) + theme(legend.text = element_text(size=22)) + theme(axis.title.x = element_text(size=26,colour = "black"),axis.text.x=element_text(angle=270, vjust=0.5, hjust=0.5, size=22))+theme(axis.title.y = element_text(size=26,colour = "black"),axis.text.y=element_text(size=22)) + guides(fill = guide_legend(keywidth = 3, keyheight = 3)) 

cplot <- cplot + theme(legend.key.height=unit(5,"line")) + theme(legend.key.width=unit(3,"line")) + theme (plot.margin = unit(c(1.5,1.0,1.5,1.5), "cm")) + theme(plot.title=element_text(hjust=0.5, vjust=2.5))+ theme(axis.title.x=element_text(hjust=0.5, vjust=-1.2))+ theme(axis.title.y=element_text(hjust=0.5, vjust=0.0)) + theme(legend.title = element_text(face="plain", size = 26))

assign(paste("pdata_finalplot", d, sep=""),cplot)
print(cplot)
assign(paste("pdata_finalplotloop",d, sep=""),recordPlot())
pdata_finalplotList <- list(pdata_finalplotList, get(paste("pdata_finalplotloop", d, sep="")))
graphics.off()
}

pdf(file = paste("Missing Variables/Variables - Domains/Clustering with Missing Variables - ", MinCN, " to ", MaxCN," Clusters - Domains.pdf", sep=""), title = paste("Clustering with Missing Variables Domains", sep=""), family='Courier', width=30, height=20)
pdata_finalplotList
graphics.off()

pdata_mfinalplotList <- list()

for (y in 1:UGroups)
{
toplot<-get(paste("pdata_subset",y, sep=""))
toplot_m<-subset(toplot,variable =='Mean')
MinAxis <- min(toplot_m$Missing_Variable)
MaxAxis <- max(toplot_m$Missing_Variable)
AxisLength <- (MaxAxis-MinAxis)+1
newvalues <-c(1:AxisLength)
curvalues <- c(MinAxis:MaxAxis)
toplot_m$newVar <- newvalues[match(toplot_m$Missing_Variable, curvalues)]
toplothead <-unique(toplot_m$group)

cplot_m <-assign(paste("pdata_m_plot",y, sep=""), ggplot(toplot_m, aes(variable, value, colour = group)) +scale_x_discrete(name="Missing Variable Number", limits=c(1:AxisLength),labels=c(MinAxis:MaxAxis)) + scale_y_continuous(name="Within Cluster Sum of Squares", limits=c(MMinValue, MMaxValue), breaks = seq(MMinValue, MMaxValue, MMaxFinal)))

cplot_m <- cplot_m +geom_line(data=subset(toplot_m,variable =='Mean'), alpha=1, aes(newVar, value, size= ''), colour="gray20", inherit.aes = FALSE)+ guides (size = guide_legend(title = "Mean", override.aes = list(size = 3, colour = "gray20")))

cplot_m <- cplot_m + geom_point(data=subset(toplot_m,variable =='Mean'), aes(newVar, value), colour = "gray20", size = 10.0, alpha=1.0, inherit.aes = FALSE) + labs (shape = "Clusters") + guides(scale_alpha(guide='none')) + guides(shape = guide_legend(override.aes = list(size = 10)))

cplot_m <- cplot_m +geom_line(colour="grey20", data=subset(toplot_m,variable =='Mean'),size = 3 , alpha=1, aes(newVar, value), inherit.aes = FALSE)

cplot_m <- cplot_m + ggtitle(paste("Missing Variables Mean WCSS Values -", toplothead, "Domain", sep=" ")) + theme(plot.title = element_text(size = 40, colour = "black", face = "bold"))+ theme(legend.title = element_text(size=26)) + theme(legend.text = element_text(size=22)) + theme(axis.title.x = element_text(size=26,colour = "black"),axis.text.x=element_text(angle=270, vjust=0.5, hjust=0.5, size=22))+theme(axis.title.y = element_text(size=26,colour = "black"),axis.text.y=element_text(size=22)) + guides(fill = guide_legend(keywidth = 3, keyheight = 3)) 

cplot_m <- cplot_m + theme(legend.key.height=unit(5,"line")) + theme(legend.key.width=unit(3,"line")) + theme (plot.margin = unit(c(1.5,1.0,1.5,1.5), "cm")) + theme(plot.title=element_text(hjust=0.5, vjust=2.5))+ theme(axis.title.x=element_text(hjust=0.5, vjust=-1.2))+ theme(axis.title.y=element_text(hjust=0.5, vjust=0.0)) + theme(legend.title = element_text(face="plain", size = 26)) + guides(fill=guide_legend(title="Variable Domain"))+ theme(legend.position="none")

assign(paste("pdata_mfinalplot", y, sep=""),cplot_m)
print(cplot_m)
assign(paste("pdata_mfinalplotloop",y, sep=""),recordPlot())
pdata_mfinalplotList <- list(pdata_mfinalplotList, get(paste("pdata_mfinalplotloop", y, sep="")))
graphics.off()
}

pdf(file = paste("Missing Variables/Variables - Domains/Missing Variables Mean WCSS Values - ", MinCN, " to ", MaxCN," Clusters - Domains.pdf", sep=""), title = paste("Missing Variables Mean WCSS Values - Domains", sep=""), family='Courier', width=30, height=20)
pdata_mfinalplotList
graphics.off()

PlotType1 <- ggplot(VFinalAll, aes(Variable, Value, colour = Group)) + geom_polygon(colour=NA, aes(fill=Group), alpha=0.5) +scale_x_discrete(name="Missing Variable Number", limits=c(1:NumV)) + scale_y_continuous(name="Within Cluster Sum of Squares", limits=c(AMinValue, AMaxValue), breaks = seq(AMinValue, AMaxValue, AMaxFinal)) + scale_fill_manual(values = c('Demographic' = C1, 'Household Composition' = C2, 'Housing' = C3, 'Socio-Economic' = C4, 'Employment' = C5), breaks = c("Demographic", "Household Composition", "Housing", "Socio-Economic", "Employment"))

PlotType1 <- PlotType1 +geom_line(data=subset(pdata,variable =='Mean'), alpha=1, aes(Missing_Variable, value, size= ''), colour="gray20", inherit.aes = FALSE)+ guides (size = guide_legend(title = "Mean", override.aes = list(size = 3, colour = "gray20")))

PlotType1 <- PlotType1 +geom_line(colour="grey20", data=subset(pdata,variable =='Mean'),size = 3 , alpha=1, aes(Missing_Variable, value), inherit.aes = FALSE)

PlotType1 <- PlotType1 +geom_line(data=subset(pdata,variable!='Mean'),size=0.5, alpha=0.5, aes(Missing_Variable, value, shape=variable), colour='black', inherit.aes = FALSE) 

PlotType1 <- PlotType1 +geom_point(data=subset(pdata,variable!='Mean'), aes(Missing_Variable, value, shape = variable, col=group), size = 5.0, alpha=1.0, inherit.aes = FALSE) + labs (shape = "Clusters") + guides(scale_alpha(guide='none')) + guides(colour = "none") + guides(shape = guide_legend(override.aes = list(size = 10))) + scale_colour_manual(values = c('Demographic' = C1, 'Household Composition' = C2, 'Housing' = C3, 'Socio-Economic' = C4, 'Employment' = C5))  + theme(legend.key.height=unit(6,"line"))

PlotType1 <- PlotType1 + ggtitle("Clusters with Missing Variables") + theme(plot.title = element_text(size = 40, colour = "black", face = "bold"))+ theme(legend.title = element_text(size=26)) + theme(legend.text = element_text(size=22)) + theme(axis.title.x = element_text(size=26,colour = "black"),axis.text.x=element_text(angle=0, vjust=0.0, hjust=0.5, size=22))+theme(axis.title.y = element_text(size=26,colour = "black"),axis.text.y=element_text(size=22)) + guides(fill = guide_legend(keywidth = 3, keyheight = 3)) 

PlotType1 <- PlotType1 + theme(legend.key.height=unit(5,"line")) + theme(legend.key.width=unit(3,"line")) + theme (plot.margin = unit(c(1.5,1.0,1.5,1.5), "cm")) + theme(plot.title=element_text(hjust=0.5, vjust=2.5))+ theme(axis.title.x=element_text(hjust=0.5, vjust=-1.2))+ theme(axis.title.y=element_text(hjust=0.5, vjust=0.0)) + theme(legend.title = element_text(face="plain", size = 26)) + guides(fill=guide_legend(title="Variable Domain"))

PlotType2 <- ggplot(VFinalAll, aes(Variable, Value, colour = Group)) + geom_polygon(colour=NA, aes(fill=Group), alpha=0.5) +scale_x_discrete(name="Missing Variable Number", limits=c(1:NumV)) + scale_y_continuous(name="Within Cluster Sum of Squares", limits=c(AMinValue, AMaxValue), breaks = seq(AMinValue, AMaxValue, AMaxFinal)) + scale_fill_manual(values = c('Demographic' = C1, 'Household Composition' = C2, 'Housing' = C3, 'Socio-Economic' = C4, 'Employment' = C5), breaks = c("Demographic", "Household Composition", "Housing", "Socio-Economic", "Employment"))

PlotType2 <- PlotType2 +geom_line(data=subset(pdata,variable =='Mean'), alpha=1, aes(Missing_Variable, value, size= ''), colour="gray20", inherit.aes = FALSE)+ guides (size = guide_legend(title = "Mean", override.aes = list(size = 3, colour = "gray20")))

PlotType2 <- PlotType2 +geom_line(colour="grey20", data=subset(pdata,variable =='Mean'),size = 3 , alpha=1, aes(Missing_Variable, value), inherit.aes = FALSE)

PlotType2 <- PlotType2 +geom_line(data=subset(pdata,variable!='Mean'),size=0.5, alpha=0.5, aes(Missing_Variable, value, shape=variable), colour='black', inherit.aes = FALSE) 

PlotType2 <- PlotType2 +geom_point(data=subset(pdata,variable!='Mean'), aes(Missing_Variable, value, shape = variable, col=group), size = 5.0, alpha=1.0, inherit.aes = FALSE) + labs (shape = "Clusters") + guides(scale_alpha(guide='none')) + guides(colour = "none") + guides(shape = guide_legend(override.aes = list(size = 10))) + scale_colour_manual(values = c('Demographic' = C1, 'Household Composition' = C2, 'Housing' = C3, 'Socio-Economic' = C4, 'Employment' = C5))  + theme(legend.key.height=unit(6,"line"))

PlotType2 <- PlotType2 + ggtitle("Clusters with Missing Variables") + theme(plot.title = element_text(size = 40, colour = "black", face = "bold"))+ theme(legend.title = element_text(size=26)) + theme(legend.text = element_text(size=22)) + theme(axis.title.x = element_text(size=26,colour = "black"),axis.text.x=element_text(angle=270, vjust=0.5, hjust=0.5, size=22))+theme(axis.title.y = element_text(size=26,colour = "black"),axis.text.y=element_text(size=22)) + guides(fill = guide_legend(keywidth = 3, keyheight = 3)) 

PlotType2 <- PlotType2 + theme(legend.key.height=unit(5,"line")) + theme(legend.key.width=unit(3,"line")) + theme (plot.margin = unit(c(1.5,1.0,1.5,1.5), "cm")) + theme(plot.title=element_text(hjust=0.5, vjust=2.5))+ theme(axis.title.x=element_text(hjust=0.5, vjust=-1.2))+ theme(axis.title.y=element_text(hjust=0.5, vjust=0.0)) + theme(legend.title = element_text(face="plain", size = 26)) + guides(fill=guide_legend(title="Variable Domain"))

PlotType3 <- ggplot(VFinalAll, aes(Variable, Value, colour = Group)) + geom_polygon(colour=NA, aes(fill=Group), alpha=0.5) +scale_x_discrete(name="Missing Variable", limits=c(1:NumV), breaks = seq(1, NumV, 1), labels =c(1, rep("",NumV-2), NumV)) + scale_y_continuous(name="Within Cluster Sum of Squares", limits=c(AMinValue, AMaxValue), breaks = seq(AMinValue, AMaxValue, AMaxFinal)) + scale_fill_manual(values = c('Demographic' = C1, 'Household Composition' = C2, 'Housing' = C3, 'Socio-Economic' = C4, 'Employment' = C5), breaks = c("Demographic", "Household Composition", "Housing", "Socio-Economic", "Employment")) + theme(axis.text.x=element_blank())

PlotType3 <- PlotType3 +geom_line(data=subset(pdata,variable =='Mean'), alpha=1, aes(Missing_Variable, value, size= ''), colour="gray20", inherit.aes = FALSE)+ guides (size = guide_legend(title = "Mean", override.aes = list(size = 3, colour = "gray20")))

PlotType3 <- PlotType3 +geom_line(colour="grey20", data=subset(pdata,variable =='Mean'),size = 3 , alpha=1, aes(Missing_Variable, value), inherit.aes = FALSE)

PlotType3 <- PlotType3 +geom_line(data=subset(pdata,variable!='Mean'),size=1.5, alpha=0.5, aes(Missing_Variable, value, linetype=variable), colour='black', inherit.aes = FALSE)+scale_linetype_discrete(name = "Clusters")

PlotType3 <- PlotType3 +ggtitle("Clusters with Missing Variables") + theme(plot.title = element_text(size = 40, colour = "black", face = "bold"))+ theme(legend.title = element_text(size=26)) + theme(legend.text = element_text(size=22)) + theme(axis.title.x = element_text(size=26,colour = "black"))+theme(axis.title.y = element_text(size=26,colour = "black"),axis.text.x=element_text(size=22),axis.text.y=element_text(size=22)) + guides(fill = guide_legend(keywidth = 3, keyheight = 3)) 

PlotType3 <- PlotType3 + theme(legend.key.height=unit(5,"line")) + theme(legend.key.width=unit(3,"line")) + theme (plot.margin = unit(c(1.5,1.0,1.5,1.5), "cm")) + theme(plot.title=element_text(hjust=0.5, vjust=2.5))+ theme(axis.title.x=element_text(hjust=0.5, vjust=-1.2))+ theme(axis.title.y=element_text(hjust=0.5, vjust=0.0)) + theme(legend.title = element_text(face="plain", size = 26)) + guides(fill=guide_legend(title="Variable Domain"))

PlotType4 <- ggplot(VFinalAll, aes(Variable, Value, colour = Group)) + geom_polygon(colour=NA, aes(fill=Group), alpha=0.5) +scale_x_discrete(name="Missing Variable Number", limits=c(1:NumV)) + scale_y_continuous(name="Within Cluster Sum of Squares", limits=c(AMinValue, AMaxValue), breaks = seq(AMinValue, AMaxValue, AMaxFinal)) + scale_fill_manual(values = c('Demographic' = C1, 'Household Composition' = C2, 'Housing' = C3, 'Socio-Economic' = C4, 'Employment' = C5), breaks = c("Demographic", "Household Composition", "Housing", "Socio-Economic", "Employment"))

PlotType4 <- PlotType4 +geom_line(data=subset(pdata,variable =='Mean'), alpha=1, aes(Missing_Variable, value, size= ''), colour="gray20", inherit.aes = FALSE)+ guides (size = guide_legend(title = "Mean", override.aes = list(size = 3, colour = "gray20")))

PlotType4 <- PlotType4 +geom_line(colour="grey20", data=subset(pdata,variable =='Mean'),size = 3 , alpha=1, aes(Missing_Variable, value), inherit.aes = FALSE)

PlotType4 <- PlotType4 +geom_line(data=subset(pdata,variable!='Mean'),size=1.5, alpha=0.5, aes(Missing_Variable, value, linetype=variable), colour='black', inherit.aes = FALSE)+scale_linetype_discrete(name = "Clusters")

PlotType4 <- PlotType4 +ggtitle("Clusters with Missing Variables") + theme(plot.title = element_text(size = 40, colour = "black", face = "bold"))+ theme(legend.title = element_text(size=26)) + theme(legend.text = element_text(size=22)) + theme(axis.title.x = element_text(size=26,colour = "black"))+theme(axis.title.y = element_text(size=26,colour = "black"),axis.text.x=element_text(angle=270, vjust=0.5, hjust=0.5, size=22))+theme(axis.title.y = element_text(size=26,colour = "black"),axis.text.y=element_text(size=22)) + guides(fill = guide_legend(keywidth = 3, keyheight = 3)) 

PlotType4 <- PlotType4 + theme(legend.key.height=unit(5,"line")) + theme(legend.key.width=unit(3,"line")) + theme (plot.margin = unit(c(1.5,1.0,1.5,1.5), "cm")) + theme(plot.title=element_text(hjust=0.5, vjust=2.5))+ theme(axis.title.x=element_text(hjust=0.5, vjust=-1.2))+ theme(axis.title.y=element_text(hjust=0.5, vjust=0.0)) + theme(legend.title = element_text(face="plain", size = 26)) + guides(fill=guide_legend(title="Variable Domain"))

PlotTypeMean1 <- ggplot(VFinal_MeanAll, aes(Variable, Value, colour = Group)) + geom_polygon(colour=NA, aes(fill=Group), alpha=0.5) +scale_x_discrete(name="Missing Variable Number", limits=c(1:NumV)) + scale_y_continuous(name="Within Cluster Sum of Squares", limits=c(MMinValue, MMaxValue), breaks = seq(MMinValue, MMaxValue, MMaxFinal)) + scale_fill_manual(values = c('Demographic' = C1, 'Household Composition' = C2, 'Housing' = C3, 'Socio-Economic' = C4, 'Employment' = C5), breaks = c("Demographic", "Household Composition", "Housing", "Socio-Economic", "Employment"))

PlotTypeMean1 <- PlotTypeMean1 +geom_line(data=subset(pdata_m,variable =='Mean'), alpha=1, aes(Missing_Variable, value, size= ''), colour="gray20", inherit.aes = FALSE)+ guides (size = guide_legend(title = "Mean", override.aes = list(size = 3, colour = "gray20")))

PlotTypeMean1 <- PlotTypeMean1 + geom_point(data=subset(pdata_m,variable =='Mean'), aes(Missing_Variable, value), colour = "gray20", size = 10.0, alpha=1.0, inherit.aes = FALSE) + labs (shape = "Clusters") + guides(scale_alpha(guide='none')) + guides(shape = guide_legend(override.aes = list(size = 10)))

PlotTypeMean1 <- PlotTypeMean1 +geom_line(colour="grey20", data=subset(pdata_m,variable =='Mean'),size = 3 , alpha=1, aes(Missing_Variable, value), inherit.aes = FALSE)

PlotTypeMean1 <- PlotTypeMean1 + ggtitle("Mean Within Cluster Sum of Squares Values for Missing Variables") + theme(plot.title = element_text(size = 40, colour = "black", face = "bold"))+ theme(legend.title = element_text(size=26)) + theme(legend.text = element_text(size=22)) + theme(axis.title.x = element_text(size=26,colour = "black"),axis.text.x=element_text(angle=0, vjust=0.0, hjust=0.5, size=22))+theme(axis.title.y = element_text(size=26,colour = "black"),axis.text.y=element_text(size=22)) + guides(fill = guide_legend(keywidth = 3, keyheight = 3)) 

PlotTypeMean1 <- PlotTypeMean1 + theme(legend.key.height=unit(5,"line")) + theme(legend.key.width=unit(3,"line")) + theme (plot.margin = unit(c(1.5,1.0,1.5,1.5), "cm")) + theme(plot.title=element_text(hjust=0.5, vjust=2.5))+ theme(axis.title.x=element_text(hjust=0.5, vjust=-1.2))+ theme(axis.title.y=element_text(hjust=0.5, vjust=0.0)) + theme(legend.title = element_text(face="plain", size = 26)) + guides(fill=guide_legend(title="Variable Domain"))

PlotTypeMean2 <- ggplot(VFinal_MeanAll, aes(Variable, Value, colour = Group)) + geom_polygon(colour=NA, aes(fill=Group), alpha=0.5) +scale_x_discrete(name="Missing Variable Number", limits=c(1:NumV)) + scale_y_continuous(name="Within Cluster Sum of Squares", limits=c(MMinValue, MMaxValue), breaks = seq(MMinValue, MMaxValue, MMaxFinal)) + scale_fill_manual(values = c('Demographic' = C1, 'Household Composition' = C2, 'Housing' = C3, 'Socio-Economic' = C4, 'Employment' = C5), breaks = c("Demographic", "Household Composition", "Housing", "Socio-Economic", "Employment"))

PlotTypeMean2 <- PlotTypeMean2 +geom_line(data=subset(pdata_m,variable =='Mean'), alpha=1, aes(Missing_Variable, value, size= ''), colour="gray20", inherit.aes = FALSE)+ guides (size = guide_legend(title = "Mean", override.aes = list(size = 3, colour = "gray20")))

PlotTypeMean2 <- PlotTypeMean2 + geom_point(data=subset(pdata_m,variable =='Mean'), aes(Missing_Variable, value), colour = "gray20", size = 10.0, alpha=1.0, inherit.aes = FALSE) + labs (shape = "Clusters") + guides(scale_alpha(guide='none')) + guides(shape = guide_legend(override.aes = list(size = 10)))

PlotTypeMean2 <- PlotTypeMean2 +geom_line(colour="grey20", data=subset(pdata_m,variable =='Mean'),size = 3 , alpha=1, aes(Missing_Variable, value), inherit.aes = FALSE)

PlotTypeMean2 <- PlotTypeMean2 + ggtitle("Mean Within Cluster Sum of Squares Values for Missing Variables") + theme(plot.title = element_text(size = 40, colour = "black", face = "bold"))+ theme(legend.title = element_text(size=26)) + theme(legend.text = element_text(size=22)) + theme(axis.title.x = element_text(size=26,colour = "black"),axis.text.x=element_text(angle=270, vjust=0.5, hjust=0.5, size=22))+theme(axis.title.y = element_text(size=26,colour = "black"),axis.text.y=element_text(size=22)) + guides(fill = guide_legend(keywidth = 3, keyheight = 3)) 

PlotTypeMean2 <- PlotTypeMean2 + theme(legend.key.height=unit(5,"line")) + theme(legend.key.width=unit(3,"line")) + theme (plot.margin = unit(c(1.5,1.0,1.5,1.5), "cm")) + theme(plot.title=element_text(hjust=0.5, vjust=2.5))+ theme(axis.title.x=element_text(hjust=0.5, vjust=-1.2))+ theme(axis.title.y=element_text(hjust=0.5, vjust=0.0)) + theme(legend.title = element_text(face="plain", size = 26)) + guides(fill=guide_legend(title="Variable Domain"))

PlotTypeMean3 <- ggplot(VFinal_MeanAll, aes(Variable, Value, colour = Group)) + geom_polygon(colour=NA, aes(fill=Group), alpha=0.5) +scale_x_discrete(name="Missing Variable", limits=c(1:NumV), breaks = seq(1, NumV, 1), labels =c(1, rep("",NumV-2), NumV)) + scale_y_continuous(name="Within Cluster Sum of Squares", limits=c(MMinValue, MMaxValue), breaks = seq(MMinValue, MMaxValue, MMaxFinal)) + scale_fill_manual(values = c('Demographic' = C1, 'Household Composition' = C2, 'Housing' = C3, 'Socio-Economic' = C4, 'Employment' = C5), breaks = c("Demographic", "Household Composition", "Housing", "Socio-Economic", "Employment")) + theme(axis.text.x=element_blank())

PlotTypeMean3 <- PlotTypeMean3 +geom_line(data=subset(pdata_m,variable =='Mean'), alpha=1, aes(Missing_Variable, value, size= ''), colour="gray20", inherit.aes = FALSE)+ guides (size = guide_legend(title = "Mean", override.aes = list(size = 3, colour = "gray20")))

PlotTypeMean3 <- PlotTypeMean3 + geom_point(data=subset(pdata_m,variable =='Mean'), aes(Missing_Variable, value), colour = "gray20", size = 10.0, alpha=0.0, inherit.aes = FALSE) + labs (shape = "Clusters") + guides(scale_alpha(guide='none')) + guides(shape = guide_legend(override.aes = list(size = 10)))

PlotTypeMean3 <- PlotTypeMean3 +geom_line(colour="grey20", data=subset(pdata_m,variable =='Mean'),size = 3 , alpha=1, aes(Missing_Variable, value), inherit.aes = FALSE)

PlotTypeMean3 <- PlotTypeMean3 + ggtitle("Mean Within Cluster Sum of Squares Values for Missing Variables")+ theme(plot.title = element_text(size = 40, colour = "black", face = "bold"))+ theme(legend.title = element_text(size=26)) + theme(legend.text = element_text(size=22)) + theme(axis.title.x = element_text(size=26,colour = "black"),axis.text.x=element_text(size=22))+theme(axis.title.y = element_text(size=26,colour = "black"),axis.text.y=element_text(size=22)) + guides(fill = guide_legend(keywidth = 3, keyheight = 3)) 

PlotTypeMean3 <- PlotTypeMean3 + theme(legend.key.height=unit(5,"line")) + theme(legend.key.width=unit(3,"line")) + theme (plot.margin = unit(c(1.5,1.0,1.5,1.5), "cm")) + theme(plot.title=element_text(hjust=0.5, vjust=2.5))+ theme(axis.title.x=element_text(hjust=0.5, vjust=-1.2))+ theme(axis.title.y=element_text(hjust=0.5, vjust=0.0)) + theme(legend.title = element_text(face="plain", size = 26)) + guides(fill=guide_legend(title="Variable Domain"))

PlotTypeMean4 <- ggplot(VFinal_MeanAll, aes(Variable, Value, colour = Group)) + geom_polygon(colour=NA, aes(fill=Group), alpha=0.5) +scale_x_discrete(name="Missing Variable Number", limits=c(1:NumV)) + scale_y_continuous(name="Within Cluster Sum of Squares", limits=c(MMinValue, MMaxValue), breaks = seq(MMinValue, MMaxValue, MMaxFinal)) + scale_fill_manual(values = c('Demographic' = C1, 'Household Composition' = C2, 'Housing' = C3, 'Socio-Economic' = C4, 'Employment' = C5), breaks = c("Demographic", "Household Composition", "Housing", "Socio-Economic", "Employment"))

PlotTypeMean4 <- PlotTypeMean4 +geom_line(data=subset(pdata_m,variable =='Mean'), alpha=1, aes(Missing_Variable, value, size= ''), colour="gray20", inherit.aes = FALSE)+ guides (size = guide_legend(title = "Mean", override.aes = list(size = 3, colour = "gray20")))

PlotTypeMean4 <- PlotTypeMean4 + geom_point(data=subset(pdata_m,variable =='Mean'), aes(Missing_Variable, value), colour = "gray20", size = 10.0, alpha=0.0, inherit.aes = FALSE) + labs (shape = "Clusters") + guides(scale_alpha(guide='none')) + guides(shape = guide_legend(override.aes = list(size = 10)))

PlotTypeMean4 <- PlotTypeMean4 +geom_line(colour="grey20", data=subset(pdata_m,variable =='Mean'),size = 3 , alpha=1, aes(Missing_Variable, value), inherit.aes = FALSE)

PlotTypeMean4 <- PlotTypeMean4 + ggtitle("Mean Within Cluster Sum of Squares Values for Missing Variables")+ theme(plot.title = element_text(size = 40, colour = "black", face = "bold"))+ theme(legend.title = element_text(size=26)) + theme(legend.text = element_text(size=22)) + theme(axis.title.x = element_text(size=26,colour = "black"),axis.text.x=element_text(angle=270, vjust=0.5, hjust=0.5, size=22))+theme(axis.title.y = element_text(size=26,colour = "black"),axis.text.y=element_text(size=22)) + guides(fill = guide_legend(keywidth = 3, keyheight = 3)) 

PlotTypeMean4 <- PlotTypeMean4 + theme(legend.key.height=unit(5,"line")) + theme(legend.key.width=unit(3,"line")) + theme (plot.margin = unit(c(1.5,1.0,1.5,1.5), "cm")) + theme(plot.title=element_text(hjust=0.5, vjust=2.5))+ theme(axis.title.x=element_text(hjust=0.5, vjust=-1.2))+ theme(axis.title.y=element_text(hjust=0.5, vjust=0.0)) + theme(legend.title = element_text(face="plain", size = 26)) + guides(fill=guide_legend(title="Variable Domain"))

PlotType_List<- list(PlotType1, PlotType2, PlotType3, PlotType4)
pdf(file = paste("Missing Variables/Variables - All/Clustering with Missing Variables - ", MinCN, " to ", MaxCN," Clusters - 4 Plot Types.pdf", sep=""), title = paste("Clustering with Missing Variables - ", MinCN, " to ", MaxCN," Clusters - 4 Plot Types"), family='Courier', width=30, height=20)
PlotType_List
graphics.off()

PlotTypeMean_List<- list(PlotTypeMean1, PlotTypeMean2, PlotTypeMean3, PlotTypeMean4)
pdf(file = paste("Missing Variables/Variables - All/Missing Variables Mean WCSS Values - ", MinCN, " to ", MaxCN," Clusters - 4 Plot Types.pdf", sep=""), title = paste("Missing Variables Mean WCSS Values - ", MinCN, " to ", MaxCN," Clusters - 4 Plot Types"), family='Courier', width=30, height=20)
PlotTypeMean_List
graphics.off()

if (RQCLUSUP == "YES")
{
unlink("Missing Variables/Cluster Update", recursive = TRUE,force = TRUE)
}

MissingEnd <- Sys.time()

CorrelationEnd-CorrelationStart
MissingEnd-MissingStart
MissingEnd-CorrelationStart
EndPlots-StartPlots

##########SAVE WORKSPACE##########
gc()
save.image(paste("2011 OAC - Selecting Variables - ",ColNumber, " - All Correlation and Missing Variable Clustering Data and Outputs.RData", sep=""))
##################################

