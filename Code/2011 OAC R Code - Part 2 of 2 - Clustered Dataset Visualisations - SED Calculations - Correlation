#####    #####     #####    #####     #####    #####     #####    #####     #####    #####     #####
# Code to Create the 2011 Area Classification for Output Areas (2011 OAC)
# Feel free to share and reuse with attribution
# Part 2 of 2 - Clustered Dataset Visualisations, SED Calculations and Correlation
#####    #####     #####    #####     #####    #####     #####    #####     #####    #####     #####

####################################################################################################
# Setup ############################################################################################
####################################################################################################

#Library packages required
library(modeest)
library(reshape2)
library(plyr)
library(ggplot2)
library(RColorBrewer)
library(plotrix)
library(fpc)
library(ineq)
library(grid)
library(Hmisc)
library(car)
library(gtools)

#Set Working Directory
#setwd("/Users/Chris/Documents/2011 OAC")

#Do you wish to have the SED values calculated?
#Enter "YES"/"NO"
RQSED <- "YES"

#Do you wish to use 2011 OAC titles and document names?
#Enter "YES"/"NO"
RQOAC <- "YES"

#Number of columns for plots
FCol <- 9

#Number of rows for plots
FRow <- 7

OutputsStart <- Sys.time()

####################################################################################################
# Bar Plots Setup ##################################################################################
####################################################################################################

if(RQCLUS=="SUPERGROUP")
{
BARCOLOURHIGH<-"midnightblue"
BARCOLOURLOW<-"skyblue3"
RADIALCOLOUR<-"blue"
}

if(RQCLUS=="GROUP")
{
BARCOLOURHIGH<-"darkgreen"
BARCOLOURLOW<-"darkseagreen3"
RADIALCOLOUR<-"forestgreen"
}

if(RQCLUS=="SUBGROUP")
{
BARCOLOURHIGH<-"purple4"
BARCOLOURLOW<-"plum3"
RADIALCOLOUR<-"purple4"
}

BARPLOT_B <- "Cluster Profiles"
BARPLOT_V <- "Cluster"
RADIALNAME_B <- "Cluster Profiles"
RADIALSUB_B <- "Cluster Profile"
DOCFILE_B <- "Cluster"
BARSUB_G <- "Global Mean"
BARSUB_C <- "Parent Cluster Mean"
BARSUB_CG <- "Global and Parent Cluster Mean"

if(RQOAC =="YES")
{
	if(RQCLUS=="SUPERGROUP")
	{
	BARPLOT_B <- "2011 OAC Supergroups"
	BARPLOT_V <- "Supergroup"
	RADIALNAME_B <- "2011 OAC Supergroups"
	RADIALSUB_B <- "2011 OAC Supergroup Profile"
	DOCFILE_B <- "2011 OAC Supergroup"
	}
	if(RQCLUS=="GROUP")
	{
	BARPLOT_B <- "2011 OAC Groups"
	BARPLOT_V <- "Group"
	RADIALNAME_B <- "2011 OAC Groups"
	RADIALSUB_B <- "2011 OAC Group Profile"
	DOCFILE_B <- "2011 OAC Group"
	}
	if(RQCLUS=="SUBGROUP")
	{
	BARPLOT_B <- "2011 OAC Subgroups"
	BARPLOT_V <- "Subgroup"
	RADIALNAME_B <- "2011 OAC Subgroups"
	RADIALSUB_B <- "2011 OAC Subgroup Profile"
	DOCFILE_B <- "2011 OAC Subgroup"
	}
}

####################################################################################################
# Bar Plots Parent Cluster Mean ####################################################################
####################################################################################################

BARPLOTPARENTHIGH = BARCOLOURHIGH
BARPLOTPARENTLOW = BARCOLOURLOW
if(RQCLUS!="SUPERGROUP")
{
BARPLOTPARENTHIGH = BARCOLOURLOW
BARPLOTPARENTLOW = BARCOLOURHIGH
}

dir.create("Cluster Plots", showWarnings = FALSE)

OAC_Converted_Transformed_Range_Clus_Nu<- nrow(data.frame(unique(OAC_Converted_Transformed_Range$Cluster)))

OAC_Converted_Transformed_Range_Global_Mean<- data.frame(lapply(OAC_Converted_Transformed_Range[,1:c(ncol(OAC_Converted_Transformed_Range)-1)], mean))

OAC_Converted_Transformed_Range_Out_Tab<- data.frame(names(OAC_Converted_Transformed_Range_Global_Mean))
colnames(OAC_Converted_Transformed_Range_Out_Tab)<- c("OAC_Converted_Transformed_Range_Var")

for (i in 1:OAC_Converted_Transformed_Range_Clus_Nu)
	{
	OAC_Converted_Transformed_Range_Sel<- OAC_Converted_Transformed_Range[which(OAC_Converted_Transformed_Range$Cluster==i),1:c(ncol(OAC_Converted_Transformed_Range)-1)]
	if(RQCLUS=="GROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, letters[i], sep="")
	OAC_Converted_Transformed_Range_Sel<- OAC_Converted_Transformed_Range[which(OAC_Converted_Transformed_Range$Cluster==BarPlotClus),1:c(ncol(OAC_Converted_Transformed_Range)-1)]
	}
	if(RQCLUS=="SUBGROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, i, sep="")
	OAC_Converted_Transformed_Range_Sel<- OAC_Converted_Transformed_Range[which(OAC_Converted_Transformed_Range$Cluster==BarPlotClus),1:c(ncol(OAC_Converted_Transformed_Range)-1)]
	}
	
	OAC_Converted_Transformed_Range_Mean<- data.frame(lapply(OAC_Converted_Transformed_Range_Sel, mean))
	OAC_Converted_Transformed_Range_Out_Col<- t(OAC_Converted_Transformed_Range_Mean-OAC_Converted_Transformed_Range_Global_Mean)
	colnames(OAC_Converted_Transformed_Range_Out_Col)<- c(paste(BARPLOT_V, i, sep=" "))
	
	if(RQCLUS=="GROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, letters[i], sep="")
	colnames(OAC_Converted_Transformed_Range_Out_Col)<- c(paste(BARPLOT_V, BarPlotClus, sep=" "))
	}
	if(RQCLUS=="SUBGROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, i, sep="")
	colnames(OAC_Converted_Transformed_Range_Out_Col)<- c(paste(BARPLOT_V, BarPlotClus, sep=" "))
	}
	OAC_Converted_Transformed_Range_Out_Tab<- data.frame(OAC_Converted_Transformed_Range_Out_Tab, OAC_Converted_Transformed_Range_Out_Col)
	}
	
OAC_Converted_Transformed_Range_Out_Var_ID <- data.frame (1:nrow(OAC_Converted_Transformed_Range_Out_Tab),OAC_Converted_Transformed_Range_Out_Tab)
colnames(OAC_Converted_Transformed_Range_Out_Var_ID) <- c("OAC_Converted_Transformed_Range_VarID",names(OAC_Converted_Transformed_Range_Out_Tab))
colnames(OAC_Converted_Transformed_Range_Out_Var_ID) <- gsub('.', ' ',colnames(OAC_Converted_Transformed_Range_Out_Var_ID),fixed=TRUE)

OAC_Converted_Transformed_Range_Final <-NULL

for (i in 1:OAC_Converted_Transformed_Range_Clus_Nu)
	{
	if(RQCLUS=="SUPERGROUP")
	{
	OAC_Converted_Transformed_Range_LDat<- data.frame(i, OAC_Converted_Transformed_Range_Out_Var_ID$OAC_Converted_Transformed_Range_VarID, OAC_Converted_Transformed_Range_Out_Var_ID[paste(BARPLOT_V," ",i, sep="")])
	}
	if(RQCLUS=="GROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, letters[i], sep="")
	OAC_Converted_Transformed_Range_LDat<- data.frame(BarPlotClus, OAC_Converted_Transformed_Range_Out_Var_ID$OAC_Converted_Transformed_Range_VarID, OAC_Converted_Transformed_Range_Out_Var_ID[paste(BARPLOT_V," ",BarPlotClus, sep="")])
	}
	if(RQCLUS=="SUBGROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, i, sep="")
	OAC_Converted_Transformed_Range_LDat<- data.frame(BarPlotClus, OAC_Converted_Transformed_Range_Out_Var_ID$OAC_Converted_Transformed_Range_VarID, OAC_Converted_Transformed_Range_Out_Var_ID[paste(BARPLOT_V," ",BarPlotClus, sep="")])
	}
	names(OAC_Converted_Transformed_Range_LDat)<- c("clusternumber", "ID","data")
	OAC_Converted_Transformed_Range_LDat$clusterrank <- ave(OAC_Converted_Transformed_Range_LDat$data, FUN=rank)
	OAC_Converted_Transformed_Range_Final <- rbind(OAC_Converted_Transformed_Range_Final, OAC_Converted_Transformed_Range_LDat) 
	}

names(OAC_Converted_Transformed_Range_Final)<- c("clusterNumber", "varID", "Data", "clusterRank")

OAC_Converted_Transformed_Range_CPD<-dcast(OAC_Converted_Transformed_Range_Final, varID ~ clusterNumber, value.var = "Data")
colnames(OAC_Converted_Transformed_Range_CPD)[1] <-"Variable Number"
for (ucl in 1:CN)
{
	if(RQCLUS=="SUPERGROUP")
	{
	colnames(OAC_Converted_Transformed_Range_CPD)[ucl+1] <-paste(BARPLOT_V, ucl, "Variable Data", sep=" ")
	}
	if(RQCLUS=="GROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, letters[ucl], sep="")
	colnames(OAC_Converted_Transformed_Range_CPD)[ucl+1] <-paste(BARPLOT_V, BarPlotClus, "Variable Data", sep=" ")
	}
	if(RQCLUS=="SUBGROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, ucl, sep="")
	colnames(OAC_Converted_Transformed_Range_CPD)[ucl+1] <-paste(BARPLOT_V, BarPlotClus, "Variable Data", sep=" ")
	}
}

for (lse in 1:CN+1)
{

OAC_Converted_Transformed_Range_VarIDList<-OAC_Converted_Transformed_Range_CPD[1]

OAC_Converted_Transformed_RangeClusterList<-OAC_Converted_Transformed_Range_CPD[lse]
colnames(OAC_Converted_Transformed_Range_VarIDList) <- paste(BARPLOT_V, lse-1, "Variable Number", sep=" ")
	

	if(RQCLUS=="GROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, letters[lse-1], sep="")
	colnames(OAC_Converted_Transformed_Range_VarIDList) <- paste(BARPLOT_V, BarPlotClus, "Variable Number", sep=" ")
	}

	if(RQCLUS=="SUBGROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, lse-1, sep="")
	colnames(OAC_Converted_Transformed_Range_VarIDList) <- paste(BARPLOT_V, BarPlotClus, "Variable Number", sep=" ")
	}

OAC_Converted_Transformed_RangeClusterData <- cbind(OAC_Converted_Transformed_Range_VarIDList,OAC_Converted_Transformed_RangeClusterList)
OAC_Converted_Transformed_RangeClusterDataOrder<-OAC_Converted_Transformed_RangeClusterData[order(-OAC_Converted_Transformed_RangeClusterData[,2]),]
assign(paste("OAC_Converted_Transformed_RangeClusterDataOrder", lse-1, sep=""), OAC_Converted_Transformed_RangeClusterDataOrder)
}

OAC_Converted_Transformed_RangeClusterDataOrderAllCMean<-get(paste("OAC_Converted_Transformed_RangeClusterDataOrder", 1, sep=""))

for (kcl in 1:CN)
{
OAC_Converted_Transformed_RangeClusterDataOrderKCL<-get(paste("OAC_Converted_Transformed_RangeClusterDataOrder", kcl, sep=""))
HJ<- paste(BARPLOT_V, kcl, sep=" ")
	if(RQCLUS=="GROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, letters[kcl], sep="")
	HJ<- paste(BARPLOT_V, BarPlotClus, sep=" ")
	}

	if(RQCLUS=="SUBGROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, kcl, sep="")
	HJ<- paste(BARPLOT_V, BarPlotClus, sep=" ")
	}
colnames(OAC_Converted_Transformed_RangeClusterDataOrderKCL)[1]<- paste(HJ, " Variable Number",sep="")
colnames(OAC_Converted_Transformed_RangeClusterDataOrderKCL)[2]<- paste(HJ, " Variable Data",sep="")
OAC_Converted_Transformed_RangeClusterDataOrderAllCMean<-cbind(OAC_Converted_Transformed_RangeClusterDataOrderAllCMean, OAC_Converted_Transformed_RangeClusterDataOrderKCL)
}

OAC_Converted_Transformed_RangeClusterDataOrderAllCMean <- OAC_Converted_Transformed_RangeClusterDataOrderAllCMean[ ,(3:ncol(OAC_Converted_Transformed_RangeClusterDataOrderAllCMean))]

if(RQOAC =="NO")
{
write.table(OAC_Converted_Transformed_RangeClusterDataOrderAllCMean, paste("Cluster Plots/",DOCFILE_B," Profiles - Cluster Plot Data (Parent Cluster Mean).csv", sep=""), sep = ",", row.names=FALSE, col.names = TRUE, qmethod = "double")
}

if(RQOAC =="YES")
{
if(RQCLUS=="SUPERGROUP")
{
write.table(OAC_Converted_Transformed_RangeClusterDataOrderAllCMean, paste("Cluster Plots/",DOCFILE_B," Profiles - Cluster Plot Data (Global Mean).csv", sep=""), sep = ",", row.names=FALSE, col.names = TRUE, qmethod = "double")
}
if(RQCLUS=="GROUP" | RQCLUS=="SUBGROUP")
{
write.table(OAC_Converted_Transformed_RangeClusterDataOrderAllCMean, paste("Cluster Plots/",DOCFILE_B," ",Best_Cluster_Input, " Profiles - Cluster Plot Data (Parent Cluster Mean).csv", sep=""), sep = ",", row.names=FALSE, col.names = TRUE, qmethod = "double")
}
}

Data_Max<-max(OAC_Converted_Transformed_Range_Final$Data)
Data_Min<-min(OAC_Converted_Transformed_Range_Final$Data)
Plot_Max_Best_C<-round_any(Data_Max, .2, f=ceiling)
Plot_Min_Best_C<-round_any(Data_Min, -.2, f=ceiling)

cNames <- c()
for (oo in 1:CN)
	{
	cNames[oo] <- paste(BARPLOT_V, oo, sep=" ")
	
	if(RQCLUS=="GROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, letters[oo], sep="")
	cNames[oo] <- paste(BARPLOT_V, BarPlotClus, sep=" ")
	}

	if(RQCLUS=="SUBGROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, oo, sep="")
	cNames[oo] <- paste(BARPLOT_V, BarPlotClus, sep=" ")
	}
	}
	
OAC_Converted_Transformed_Range_Final$clusterName<- cNames[OAC_Converted_Transformed_Range_Final$clusterNumber]

OAC_Converted_Transformed_Range_Final$clusterName <- factor(OAC_Converted_Transformed_Range_Final$clusterName, levels = cNames)

plot<-qplot(OAC_Converted_Transformed_Range_Final$varID, OAC_Converted_Transformed_Range_Final$Data, data=OAC_Converted_Transformed_Range_Final, geom="bar", stat="identity", position = "identity",  ylim = c(Plot_Min_Best_C, Plot_Max_Best_C), xlab="Variable Number", ylab="Distance from Mean") + aes( fill=OAC_Converted_Transformed_Range_Final$clusterRank) +theme(legend.position="none")

plot<-plot + xlim (1,max(OAC_Converted_Transformed_Range_Out_Var_ID$OAC_Converted_Transformed_Range_VarID) ) + ylim (Plot_Min_Best_C, Plot_Max_Best_C) + theme(axis.text.x = element_text(angle=270, hjust = 0, vjust=0.5)) + scale_x_continuous(breaks=1:max(OAC_Converted_Transformed_Range_Out_Var_ID$OAC_Converted_Transformed_Range_VarID))

plot + scale_fill_gradient (low=brewer.pal(5,"Blues"), space="Lab", guide="colourbar")

if(RQCLUS=="SUPERGROUP")
{
BARPLOT_1 <-paste(BARPLOT_B," - ",BARSUB_G,"\n", sep="")
}

if(RQCLUS=="GROUP" | RQCLUS=="SUBGROUP")
{
BARPLOT_1 <-paste(BARPLOT_B," - ",BARSUB_C,"\n", sep="")
}

plot<-plot + facet_wrap(~ clusterName, ncol = 1) + scale_fill_gradient(low=BARPLOTPARENTLOW, high=BARPLOTPARENTHIGH, space="Lab", guide="colourbar") + ggtitle(BARPLOT_1)

plot<-plot + theme(plot.title = element_text(size = 18, colour = "black", face = "bold"))

plot<-plot + theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "white")
		)
assign(paste("ClusterMeanPlot", sep=""),plot)
print(ClusterMeanPlot)
assign(paste("ClusterMeanPlot", sep=""),recordPlot())

graphics.off()

if(RQOAC =="NO")	
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," Profiles - Bar Plots.pdf", sep=""), title = paste(RADIALNAME_B), height=16, width=11, onefile=FALSE, family='Courier')
print(ClusterMeanPlot)
graphics.off()
}

if(RQOAC =="YES")
{
if(RQCLUS=="SUPERGROUP")
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," Profiles - Bar Plots.pdf", sep=""), title = paste(RADIALNAME_B), height=16, width=11, onefile=FALSE, family='Courier')
print(ClusterMeanPlot)
graphics.off()
}
if(RQCLUS=="GROUP" | RQCLUS=="SUBGROUP")
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," ",Best_Cluster_Input, " Profiles - Bar Plots.pdf", sep=""), title = paste(RADIALNAME_B), height=16, width=11, onefile=FALSE, family='Courier')
print(ClusterMeanPlot)
graphics.off()
}
}

####################################################################################################
# Bar Plots Global Mean ############################################################################
####################################################################################################
if(exists("Global_Mean_Input")=="TRUE" && RQCLUS!="SUPERGROUP")
{
dir.create("Cluster Plots", showWarnings = FALSE)

OAC_Converted_Transformed_Range_Clus_Nu<- nrow(data.frame(unique(OAC_Converted_Transformed_Range$Cluster)))

OAC_Converted_Transformed_Range_Global_Mean<- data.frame(lapply(Global_Mean_Input[,2:c(ncol(Global_Mean_Input)-1)], mean))

OAC_Converted_Transformed_Range_Out_Tab<- data.frame(names(OAC_Converted_Transformed_Range_Global_Mean))
colnames(OAC_Converted_Transformed_Range_Out_Tab)<- c("OAC_Converted_Transformed_Range_Var")

for (i in 1:OAC_Converted_Transformed_Range_Clus_Nu)
	{
	OAC_Converted_Transformed_Range_Sel<- OAC_Converted_Transformed_Range[which(OAC_Converted_Transformed_Range$Cluster==i),1:c(ncol(OAC_Converted_Transformed_Range)-1)]
	if(RQCLUS=="GROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, letters[i], sep="")
	OAC_Converted_Transformed_Range_Sel<- OAC_Converted_Transformed_Range[which(OAC_Converted_Transformed_Range$Cluster==BarPlotClus),1:c(ncol(OAC_Converted_Transformed_Range)-1)]
	}
	if(RQCLUS=="SUBGROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, i, sep="")
	OAC_Converted_Transformed_Range_Sel<- OAC_Converted_Transformed_Range[which(OAC_Converted_Transformed_Range$Cluster==BarPlotClus),1:c(ncol(OAC_Converted_Transformed_Range)-1)]
	}
	
	OAC_Converted_Transformed_Range_Mean<- data.frame(lapply(OAC_Converted_Transformed_Range_Sel, mean))
	OAC_Converted_Transformed_Range_Out_Col<- t(OAC_Converted_Transformed_Range_Mean-OAC_Converted_Transformed_Range_Global_Mean)
	colnames(OAC_Converted_Transformed_Range_Out_Col)<- c(paste(BARPLOT_V, i, sep=" "))
	
	if(RQCLUS=="GROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, letters[i], sep="")
	colnames(OAC_Converted_Transformed_Range_Out_Col)<- c(paste(BARPLOT_V, BarPlotClus, sep=" "))
	}
	if(RQCLUS=="SUBGROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, i, sep="")
	colnames(OAC_Converted_Transformed_Range_Out_Col)<- c(paste(BARPLOT_V, BarPlotClus, sep=" "))
	}
	OAC_Converted_Transformed_Range_Out_Tab<- data.frame(OAC_Converted_Transformed_Range_Out_Tab, OAC_Converted_Transformed_Range_Out_Col)
	}
	
OAC_Converted_Transformed_Range_Out_Var_ID <- data.frame (1:nrow(OAC_Converted_Transformed_Range_Out_Tab),OAC_Converted_Transformed_Range_Out_Tab)
colnames(OAC_Converted_Transformed_Range_Out_Var_ID) <- c("OAC_Converted_Transformed_Range_VarID",names(OAC_Converted_Transformed_Range_Out_Tab))
colnames(OAC_Converted_Transformed_Range_Out_Var_ID) <- gsub('.', ' ',colnames(OAC_Converted_Transformed_Range_Out_Var_ID),fixed=TRUE)

OAC_Converted_Transformed_Range_Final <-NULL

for (i in 1:OAC_Converted_Transformed_Range_Clus_Nu)
	{
	if(RQCLUS=="SUPERGROUP")
	{
	OAC_Converted_Transformed_Range_LDat<- data.frame(i, OAC_Converted_Transformed_Range_Out_Var_ID$OAC_Converted_Transformed_Range_VarID, OAC_Converted_Transformed_Range_Out_Var_ID[paste(BARPLOT_V," ",i, sep="")])
	}
	if(RQCLUS=="GROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, letters[i], sep="")
	OAC_Converted_Transformed_Range_LDat<- data.frame(BarPlotClus, OAC_Converted_Transformed_Range_Out_Var_ID$OAC_Converted_Transformed_Range_VarID, OAC_Converted_Transformed_Range_Out_Var_ID[paste(BARPLOT_V," ",BarPlotClus, sep="")])
	}
	if(RQCLUS=="SUBGROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, i, sep="")
	OAC_Converted_Transformed_Range_LDat<- data.frame(BarPlotClus, OAC_Converted_Transformed_Range_Out_Var_ID$OAC_Converted_Transformed_Range_VarID, OAC_Converted_Transformed_Range_Out_Var_ID[paste(BARPLOT_V," ",BarPlotClus, sep="")])
	}
	names(OAC_Converted_Transformed_Range_LDat)<- c("clusternumber", "ID","data")
	OAC_Converted_Transformed_Range_LDat$clusterrank <- ave(OAC_Converted_Transformed_Range_LDat$data, FUN=rank)
	OAC_Converted_Transformed_Range_Final <- rbind(OAC_Converted_Transformed_Range_Final, OAC_Converted_Transformed_Range_LDat) 
	}

names(OAC_Converted_Transformed_Range_Final)<- c("clusterNumber", "varID", "Data", "clusterRank")

OAC_Converted_Transformed_Range_CPD<-dcast(OAC_Converted_Transformed_Range_Final, varID ~ clusterNumber, value.var = "Data")
colnames(OAC_Converted_Transformed_Range_CPD)[1] <-"Variable Number"
for (ucl in 1:CN)
{
	if(RQCLUS=="SUPERGROUP")
	{
	colnames(OAC_Converted_Transformed_Range_CPD)[ucl+1] <-paste(BARPLOT_V, ucl, "Variable Data", sep=" ")
	}
	if(RQCLUS=="GROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, letters[ucl], sep="")
	colnames(OAC_Converted_Transformed_Range_CPD)[ucl+1] <-paste(BARPLOT_V, BarPlotClus, "Variable Data", sep=" ")
	}
	if(RQCLUS=="SUBGROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, ucl, sep="")
	colnames(OAC_Converted_Transformed_Range_CPD)[ucl+1] <-paste(BARPLOT_V, BarPlotClus, "Variable Data", sep=" ")
	}
}

for (lse in 1:CN+1)
{

OAC_Converted_Transformed_Range_VarIDList<-OAC_Converted_Transformed_Range_CPD[1]

OAC_Converted_Transformed_RangeClusterList<-OAC_Converted_Transformed_Range_CPD[lse]
colnames(OAC_Converted_Transformed_Range_VarIDList) <- paste(BARPLOT_V, lse-1, "Variable Number", sep=" ")
	

	if(RQCLUS=="GROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, letters[lse-1], sep="")
	colnames(OAC_Converted_Transformed_Range_VarIDList) <- paste(BARPLOT_V, BarPlotClus, "Variable Number", sep=" ")
	}

	if(RQCLUS=="SUBGROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, lse-1, sep="")
	colnames(OAC_Converted_Transformed_Range_VarIDList) <- paste(BARPLOT_V, BarPlotClus, "Variable Number", sep=" ")
	}

OAC_Converted_Transformed_RangeClusterData <- cbind(OAC_Converted_Transformed_Range_VarIDList,OAC_Converted_Transformed_RangeClusterList)
OAC_Converted_Transformed_RangeClusterDataOrder<-OAC_Converted_Transformed_RangeClusterData[order(-OAC_Converted_Transformed_RangeClusterData[,2]),]
assign(paste("OAC_Converted_Transformed_RangeClusterDataOrder", lse-1, sep=""), OAC_Converted_Transformed_RangeClusterDataOrder)
}

OAC_Converted_Transformed_RangeClusterDataOrderAllGMean<-get(paste("OAC_Converted_Transformed_RangeClusterDataOrder", 1, sep=""))

for (kcl in 1:CN)
{
OAC_Converted_Transformed_RangeClusterDataOrderKCL<-get(paste("OAC_Converted_Transformed_RangeClusterDataOrder", kcl, sep=""))
HJ<- paste(BARPLOT_V, kcl, sep=" ")
	if(RQCLUS=="GROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, letters[kcl], sep="")
	HJ<- paste(BARPLOT_V, BarPlotClus, sep=" ")
	}

	if(RQCLUS=="SUBGROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, kcl, sep="")
	HJ<- paste(BARPLOT_V, BarPlotClus, sep=" ")
	}
colnames(OAC_Converted_Transformed_RangeClusterDataOrderKCL)[1]<- paste(HJ, " Variable Number",sep="")
colnames(OAC_Converted_Transformed_RangeClusterDataOrderKCL)[2]<- paste(HJ, " Variable Data",sep="")
OAC_Converted_Transformed_RangeClusterDataOrderAllGMean<-cbind(OAC_Converted_Transformed_RangeClusterDataOrderAllGMean, OAC_Converted_Transformed_RangeClusterDataOrderKCL)
}

OAC_Converted_Transformed_RangeClusterDataOrderAllGMean <- OAC_Converted_Transformed_RangeClusterDataOrderAllGMean[ ,(3:ncol(OAC_Converted_Transformed_RangeClusterDataOrderAllGMean))]

if(RQOAC =="NO")
{
write.table(OAC_Converted_Transformed_RangeClusterDataOrderAllGMean, paste("Cluster Plots/",DOCFILE_B," Profiles - Cluster Plot Data (Global Mean).csv", sep=""), sep = ",", row.names=FALSE, col.names = TRUE, qmethod = "double")
}

if(RQOAC =="YES")
{
if(RQCLUS=="SUPERGROUP")
{
write.table(OAC_Converted_Transformed_RangeClusterDataOrderAllGMean, paste("Cluster Plots/",DOCFILE_B," Profiles - Cluster Plot Data (Global Mean).csv", sep=""), sep = ",", row.names=FALSE, col.names = TRUE, qmethod = "double")
}
if(RQCLUS=="GROUP" | RQCLUS=="SUBGROUP")
{
write.table(OAC_Converted_Transformed_RangeClusterDataOrderAllGMean, paste("Cluster Plots/",DOCFILE_B," ",Best_Cluster_Input, " Profiles - Cluster Plot Data (Global Mean).csv", sep=""), sep = ",", row.names=FALSE, col.names = TRUE, qmethod = "double")
}
}

Data_Max<-max(OAC_Converted_Transformed_Range_Final$Data)
Data_Min<-min(OAC_Converted_Transformed_Range_Final$Data)
Plot_Max_Best_G<-round_any(Data_Max, .2, f=ceiling)
Plot_Min_Best_G<-round_any(Data_Min, -.2, f=ceiling)

cNames <- c()
for (oo in 1:CN)
	{
	cNames[oo] <- paste(BARPLOT_V, oo, sep=" ")
	
	if(RQCLUS=="GROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, letters[oo], sep="")
	cNames[oo] <- paste(BARPLOT_V, BarPlotClus, sep=" ")
	}

	if(RQCLUS=="SUBGROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, oo, sep="")
	cNames[oo] <- paste(BARPLOT_V, BarPlotClus, sep=" ")
	}
	}
	
OAC_Converted_Transformed_Range_Final$clusterName<- cNames[OAC_Converted_Transformed_Range_Final$clusterNumber]

OAC_Converted_Transformed_Range_Final$clusterName <- factor(OAC_Converted_Transformed_Range_Final$clusterName, levels = cNames)

plot<-qplot(OAC_Converted_Transformed_Range_Final$varID, OAC_Converted_Transformed_Range_Final$Data, data=OAC_Converted_Transformed_Range_Final, geom="bar", stat="identity", position = "identity",  ylim = c(Plot_Min_Best_G, Plot_Max_Best_G), xlab="Variable Number", ylab="Distance from Mean") + aes( fill=OAC_Converted_Transformed_Range_Final$clusterRank) +theme(legend.position="none")

plot<-plot + xlim (1,max(OAC_Converted_Transformed_Range_Out_Var_ID$OAC_Converted_Transformed_Range_VarID) ) + ylim (Plot_Min_Best_G, Plot_Max_Best_G) + theme(axis.text.x = element_text(angle=270, hjust = 0, vjust=0.5)) + scale_x_continuous(breaks=1:max(OAC_Converted_Transformed_Range_Out_Var_ID$OAC_Converted_Transformed_Range_VarID))

plot + scale_fill_gradient (low=brewer.pal(5,"Blues"), space="Lab", guide="colourbar")

if(RQCLUS=="SUPERGROUP")
{
BARPLOT_2 <-paste(BARPLOT_B," - ",BARSUB_G,"\n", sep="")
}

if(RQCLUS=="GROUP" | RQCLUS=="SUBGROUP")
{
BARPLOT_2 <-paste(BARPLOT_B," - ",BARSUB_G,"\n", sep="")
}

plot<-plot + facet_wrap(~ clusterName, ncol = 1) + scale_fill_gradient(low=BARCOLOURLOW, high=BARCOLOURHIGH, space="Lab", guide="colourbar") + ggtitle(BARPLOT_2)

plot<-plot + theme(plot.title = element_text(size = 18, colour = "black", face = "bold"))

plot<-plot + theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "white")
		)
assign(paste("GlobalMeanPlot", sep=""),plot)
print(GlobalMeanPlot)
assign(paste("GlobalMeanPlot", sep=""),recordPlot())
graphics.off()

if(RQOAC =="NO")	
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," Profiles - Bar Plots.pdf", sep=""), title = paste(RADIALNAME_B), height=16, width=11, onefile=TRUE, family='Courier')
}

if(RQOAC =="YES")
{
if(RQCLUS=="SUPERGROUP")
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," Profiles - Bar Plots.pdf", sep=""), title = paste(RADIALNAME_B), height=16, width=11, onefile=TRUE, family='Courier')

}
if(RQCLUS=="GROUP" | RQCLUS=="SUBGROUP")
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," ",Best_Cluster_Input, " Profiles - Bar Plots.pdf", sep=""), title = paste(RADIALNAME_B), height=16, width=11, onefile=TRUE, family='Courier')
}
}

print(GlobalMeanPlot)
print(ClusterMeanPlot)
graphics.off()
}

####################################################################################################
# Radial Plots #####################################################################################
####################################################################################################

RADIALMEANCOL <- paste("red")
if(RQCLUS!="SUPERGROUP"){RADIALMEANCOL <- paste("chocolate1")}

OAC_Converted_Transformed_Range_Cluster_Mean_Clus_Nu<- nrow(data.frame(unique(OAC_Converted_Transformed_Range$Cluster)))

OAC_Converted_Transformed_Range_Cluster_Mean_Global_Mean<- data.frame(lapply(OAC_Converted_Transformed_Range[,1:c(ncol(OAC_Converted_Transformed_Range)-1)], mean))

OAC_Converted_Transformed_Range_Cluster_Mean_Out_Tab<- data.frame(names(OAC_Converted_Transformed_Range_Cluster_Mean_Global_Mean))
colnames(OAC_Converted_Transformed_Range_Cluster_Mean_Out_Tab)<- c("OAC_Converted_Transformed_Range_Cluster_Mean_Var")

for (i in 1:OAC_Converted_Transformed_Range_Cluster_Mean_Clus_Nu)
	{
	OAC_Converted_Transformed_Range_Cluster_Mean_Sel<- OAC_Converted_Transformed_Range[which(OAC_Converted_Transformed_Range$Cluster==i),1:c(ncol(OAC_Converted_Transformed_Range)-1)]
	if(RQCLUS=="GROUP")
	{		
	BarPlotClus<-paste(Best_Cluster_Input, letters[i], sep="")
	OAC_Converted_Transformed_Range_Cluster_Mean_Sel<- OAC_Converted_Transformed_Range[which(OAC_Converted_Transformed_Range$Cluster==BarPlotClus),1:c(ncol(OAC_Converted_Transformed_Range)-1)]
	}
	if(RQCLUS=="SUBGROUP")
	{		
	BarPlotClus<-paste(Best_Cluster_Input, i, sep="")
	OAC_Converted_Transformed_Range_Cluster_Mean_Sel<- OAC_Converted_Transformed_Range[which(OAC_Converted_Transformed_Range$Cluster==BarPlotClus),1:c(ncol(OAC_Converted_Transformed_Range)-1)]
	}
	OAC_Converted_Transformed_Range_Cluster_Mean_Mean<- data.frame(lapply(OAC_Converted_Transformed_Range_Cluster_Mean_Sel, mean))
	OAC_Converted_Transformed_Range_Cluster_Mean_Out_Col<- t(OAC_Converted_Transformed_Range_Cluster_Mean_Mean-OAC_Converted_Transformed_Range_Cluster_Mean_Global_Mean)
	colnames(OAC_Converted_Transformed_Range_Cluster_Mean_Out_Col)<- c(paste("Cluster_", i, sep=""))
	if(RQCLUS=="GROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, letters[i], sep="")
	colnames(OAC_Converted_Transformed_Range_Cluster_Mean_Out_Col)<- c(paste("Cluster_", BarPlotClus, sep=""))
	}
	if(RQCLUS=="SUBGROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, i, sep="")
	colnames(OAC_Converted_Transformed_Range_Cluster_Mean_Out_Col)<- c(paste("Cluster_", BarPlotClus, sep=""))
	}
	OAC_Converted_Transformed_Range_Cluster_Mean_Out_Tab<- data.frame(OAC_Converted_Transformed_Range_Cluster_Mean_Out_Tab, OAC_Converted_Transformed_Range_Cluster_Mean_Out_Col)
	}

OAC_Converted_Transformed_Range_Cluster_Mean_Out_VarID <- data.frame (1:nrow(OAC_Converted_Transformed_Range_Cluster_Mean_Out_Tab),OAC_Converted_Transformed_Range_Cluster_Mean_Out_Tab)
colnames(OAC_Converted_Transformed_Range_Cluster_Mean_Out_VarID) <- c("OAC_Converted_Transformed_Range_Cluster_Mean_VarID",names(OAC_Converted_Transformed_Range_Cluster_Mean_Out_Tab))

OAC_Converted_Transformed_Range_Cluster_Mean_Zero_Line<- data.frame(OAC_Converted_Transformed_Range_Cluster_Mean_Out_VarID$OAC_Converted_Transformed_Range_Cluster_Mean_VarID, 0)
names(OAC_Converted_Transformed_Range_Cluster_Mean_Zero_Line)<- c("OAC_Converted_Transformed_Range_Cluster_Mean_VarID", "OAC_Converted_Transformed_Range_Cluster_Mean_Zero")

for (nn in 1:CN) {assign(paste("RPLOT", nn, sep=""), c(0))}

RPLOTList <-list()

for (ggg in 1:OAC_Converted_Transformed_Range_Cluster_Mean_Clus_Nu)
{
jjj<-ggg+2
HJ<- paste(BARPLOT_V, ggg, sep=" ")


	if(RQCLUS=="SUPERGROUP"){RADIALTITLE <- paste(HJ, sep="")}
	if(RQCLUS=="GROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, letters[ggg], sep="")
	RADIALTITLE <- paste(BARPLOT_V,BarPlotClus, sep=" ")	
	}
	if(RQCLUS=="SUBGROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, ggg, sep="")
	RADIALTITLE <- paste(BARPLOT_V,BarPlotClus, sep=" ")
	}
	
RADIALSUB_1 <-paste(RADIALSUB_B," - ",BARSUB_G, sep="")

if(RQCLUS=="GROUP" | RQCLUS=="SUBGROUP")
{
RADIALSUB_1 <-paste(RADIALSUB_B," - ",BARSUB_C, sep="")
}

par(cex.main=1.5)

radial.plot(OAC_Converted_Transformed_Range_Cluster_Mean_Zero_Line$OAC_Converted_Transformed_Range_Cluster_Mean_Zero, labels=c(OAC_Converted_Transformed_Range_Cluster_Mean_Zero_Line$OAC_Converted_Transformed_Range_Cluster_Mean_VarID), start=0,clockwise=TRUE, rp.type="p", line.col=RADIALMEANCOL,lwd=3, show.grid=TRUE, grid.col="grey", point.symbols=16,point.col="black",show.centroid=FALSE, mar=c(4,4,6,4), radial.lim=range(Plot_Min_Best_C,Plot_Max_Best_C))

radial.plot(OAC_Converted_Transformed_Range_Cluster_Mean_Out_VarID[,jjj], labels=c(OAC_Converted_Transformed_Range_Cluster_Mean_Out_VarID$OAC_Converted_Transformed_Range_Cluster_Mean_VarID), start=0,clockwise=TRUE, rp.type="p&s", line.col=RADIALCOLOUR,lwd=3, show.grid=TRUE, grid.col="grey", point.symbols=16,point.col="black",show.centroid=FALSE, radial.lim=range(Plot_Min_Best_C,Plot_Max_Best_C),  mar=c(4,4,6,4), title(main=get(paste("RADIALTITLE", sep="")), outer=FALSE, line=4), add=T)
mtext(RADIALSUB_1, line=2.5)

assign(paste("RPLOT",ggg, sep=""),recordPlot())

RPLOTList <- list(RPLOTList, get(paste("RPLOT", ggg, sep="")))

graphics.off()
}

dir.create("Cluster Plots", showWarnings = FALSE)

if(RQOAC =="NO")	
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," Profiles - Radial Plots (Individual Parent Cluster Mean).pdf", sep=""), title = RADIALNAME_B, family='Courier', width=10, height=10)
}

if(RQOAC =="YES")
{
if(RQCLUS=="SUPERGROUP")
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," Profiles - Radial Plots (Individual Global Mean).pdf", sep=""), title = RADIALNAME_B, family='Courier', width=10, height=10)
}
if(RQCLUS=="GROUP" | RQCLUS=="SUBGROUP")
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," ",Best_Cluster_Input, " Profiles - Radial Plots (Individual Parent Cluster Mean).pdf", sep=""), title = RADIALNAME_B, family='Courier', width=10, height=10)
}
}
RPLOTList
dev.off()

dir.create("Cluster Plots", showWarnings = FALSE)

if (CN<=8)
{

if(RQOAC =="NO")	
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," Profiles - Radial Plots (Combined Parent Cluster Mean).pdf", sep=""), title = RADIALNAME_B, family='Courier', width=35, height=20)
}

if(RQOAC =="YES")
{
if(RQCLUS=="SUPERGROUP")
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," Profiles - Radial Plots (Combined Global Mean).pdf", sep=""), title = RADIALNAME_B, family='Courier', width=35, height=20)
}
if(RQCLUS=="GROUP" | RQCLUS=="SUBGROUP")
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," ",Best_Cluster_Input, " Profiles - Radial Plots (Combined Parent Cluster Mean).pdf", sep=""), title = RADIALNAME_B, family='Courier', width=35, height=20)
}
}

par(mfrow = c(2,4), cex.main=2.5, xpd=TRUE, oma=c(0,0,5,0))

for (ggg in 1:OAC_Converted_Transformed_Range_Cluster_Mean_Clus_Nu)
{
jjj<-ggg+2

HJ<- paste(BARPLOT_V, ggg, sep=" ")

	if(RQCLUS=="SUPERGROUP"){RADIALTITLE <- paste(HJ, sep="")}
	if(RQCLUS=="GROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, letters[ggg], sep="")
	RADIALTITLE <- paste(BARPLOT_V,BarPlotClus, sep=" ")	
	}
	if(RQCLUS=="SUBGROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, ggg, sep="")
	RADIALTITLE <- paste(BARPLOT_V,BarPlotClus, sep=" ")
	}
	
RADIALNAME_1 <-paste(RADIALNAME_B," - ",BARSUB_G, sep="")

if(RQCLUS=="GROUP" | RQCLUS=="SUBGROUP")
{
RADIALNAME_1 <-paste(RADIALNAME_B," - ",BARSUB_C, sep="")
}

radial.plot(OAC_Converted_Transformed_Range_Cluster_Mean_Zero_Line$OAC_Converted_Transformed_Range_Cluster_Mean_Zero, labels=c(OAC_Converted_Transformed_Range_Cluster_Mean_Zero_Line$OAC_Converted_Transformed_Range_Cluster_Mean_VarID), start=0,clockwise=TRUE, rp.type="p", line.col=RADIALMEANCOL,lwd=3, show.grid=TRUE, grid.col="grey", point.symbols=16,point.col="black",show.centroid=FALSE, radial.lim=range(Plot_Min_Best_C,Plot_Max_Best_C), mar=c(5,5,15,5))

radial.plot(OAC_Converted_Transformed_Range_Cluster_Mean_Out_VarID[,jjj], labels=c(OAC_Converted_Transformed_Range_Cluster_Mean_Out_VarID$OAC_Converted_Transformed_Range_Cluster_Mean_VarID), start=0,clockwise=TRUE, rp.type="p&s", line.col=RADIALCOLOUR,lwd=3, show.grid=TRUE, grid.col="grey",
point.symbols=16,point.col="black",show.centroid=FALSE, radial.lim=range(Plot_Min_Best_C,Plot_Max_Best_C), mar=c(5,5,15,5), title(main=get(paste("RADIALTITLE", sep="")), outer=FALSE), add=T)
}
title(main=paste(RADIALNAME_1, sep=""), outer=TRUE, line=-2, cex.main=6)
graphics.off()
}

if (CN>8)

{
if(RQOAC =="NO")	
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," Profiles - Radial Plots (Combined Parent Cluster Mean).pdf", sep=""), title = RADIALNAME_B, family='Courier', width=35, height=20)
}

if(RQOAC =="YES")
{
if(RQCLUS=="SUPERGROUP")
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," Profiles - Radial Plots (Combined Global Mean).pdf", sep=""), title = RADIALNAME_B, family='Courier', width=35, height=20)
}
if(RQCLUS=="GROUP" | RQCLUS=="SUBGROUP")
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," ",Best_Cluster_Input, " Profiles - Radial Plots (Combined Parent Cluster Mean).pdf", sep=""), title = RADIALNAME_B, family='Courier', width=35, height=20)
}
}

par(mfrow = c(2,4), cex.main=2.5, xpd=FALSE, mar=c(2,2,2,2))

for (ggg in 1:OAC_Converted_Transformed_Range_Cluster_Mean_Clus_Nu)
{
jjj<-ggg+2

HJ<- paste(BARPLOT_V, ggg, sep=" ")

	if(RQCLUS=="SUPERGROUP"){RADIALTITLE <- paste(HJ, sep="")}
	if(RQCLUS=="GROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, letters[ggg], sep="")
	RADIALTITLE <- paste(BARPLOT_V,BarPlotClus, sep=" ")	
	}
	if(RQCLUS=="SUBGROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, ggg, sep="")
	RADIALTITLE <- paste(BARPLOT_V,BarPlotClus, sep=" ")
	}

radial.plot(OAC_Converted_Transformed_Range_Cluster_Mean_Zero_Line$OAC_Converted_Transformed_Range_Cluster_Mean_Zero, labels=c(OAC_Converted_Transformed_Range_Cluster_Mean_Zero_Line$OAC_Converted_Transformed_Range_Cluster_Mean_VarID), start=0,clockwise=TRUE, rp.type="p", line.col=RADIALMEANCOL,lwd=3, show.grid=TRUE, grid.col="grey", point.symbols=16,point.col="black",show.centroid=FALSE, radial.lim=range(Plot_Min_Best_C,Plot_Max_Best_C), mar=c(5,5,15,5))

radial.plot(OAC_Converted_Transformed_Range_Cluster_Mean_Out_VarID[,jjj], labels=c(OAC_Converted_Transformed_Range_Cluster_Mean_Out_VarID$OAC_Converted_Transformed_Range_Cluster_Mean_VarID), start=0,clockwise=TRUE, rp.type="p&s", line.col=RADIALCOLOUR,lwd=3, show.grid=TRUE, grid.col="grey",
point.symbols=16,point.col="black",show.centroid=FALSE, radial.lim=range(Plot_Min_Best_C,Plot_Max_Best_C), mar=c(5,5,15,5), title(main=get(paste("RADIALTITLE", sep="")), outer=FALSE), add=T)
}
graphics.off()
}

if(exists("Global_Mean_Input")=="TRUE" && RQCLUS!="SUPERGROUP")
{

OAC_Converted_Transformed_Range_Global_Mean_Clus_Nu<- nrow(data.frame(unique(OAC_Converted_Transformed_Range$Cluster)))

OAC_Converted_Transformed_Range_Global_Mean_Global_Mean<- data.frame(lapply(Global_Mean_Input[,2:c(ncol(Global_Mean_Input)-1)], mean))

OAC_Converted_Transformed_Range_Global_Mean_Out_Tab<- data.frame(names(OAC_Converted_Transformed_Range_Global_Mean_Global_Mean))
colnames(OAC_Converted_Transformed_Range_Global_Mean_Out_Tab)<- c("OAC_Converted_Transformed_Range_Global_Mean_Var")

for (i in 1:OAC_Converted_Transformed_Range_Global_Mean_Clus_Nu)
	{
	OAC_Converted_Transformed_Range_Global_Mean_Sel<- OAC_Converted_Transformed_Range[which(OAC_Converted_Transformed_Range$Cluster==i),1:c(ncol(OAC_Converted_Transformed_Range)-1)]
	if(RQCLUS=="GROUP")
	{		
	BarPlotClus<-paste(Best_Cluster_Input, letters[i], sep="")
	OAC_Converted_Transformed_Range_Global_Mean_Sel<- OAC_Converted_Transformed_Range[which(OAC_Converted_Transformed_Range$Cluster==BarPlotClus),1:c(ncol(OAC_Converted_Transformed_Range)-1)]
	}
	if(RQCLUS=="SUBGROUP")
	{		
	BarPlotClus<-paste(Best_Cluster_Input, i, sep="")
	OAC_Converted_Transformed_Range_Global_Mean_Sel<- OAC_Converted_Transformed_Range[which(OAC_Converted_Transformed_Range$Cluster==BarPlotClus),1:c(ncol(OAC_Converted_Transformed_Range)-1)]
	}
	OAC_Converted_Transformed_Range_Global_Mean_Mean<- data.frame(lapply(OAC_Converted_Transformed_Range_Global_Mean_Sel, mean))
	OAC_Converted_Transformed_Range_Global_Mean_Out_Col<- t(OAC_Converted_Transformed_Range_Global_Mean_Mean-OAC_Converted_Transformed_Range_Global_Mean_Global_Mean)
	colnames(OAC_Converted_Transformed_Range_Global_Mean_Out_Col)<- c(paste("Cluster_", i, sep=""))
	if(RQCLUS=="GROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, letters[i], sep="")
	colnames(OAC_Converted_Transformed_Range_Global_Mean_Out_Col)<- c(paste("Cluster_", BarPlotClus, sep=""))
	}
	if(RQCLUS=="SUBGROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, i, sep="")
	colnames(OAC_Converted_Transformed_Range_Global_Mean_Out_Col)<- c(paste("Cluster_", BarPlotClus, sep=""))
	}
	OAC_Converted_Transformed_Range_Global_Mean_Out_Tab<- data.frame(OAC_Converted_Transformed_Range_Global_Mean_Out_Tab, OAC_Converted_Transformed_Range_Global_Mean_Out_Col)
	}

OAC_Converted_Transformed_Range_Global_Mean_Out_VarID <- data.frame (1:nrow(OAC_Converted_Transformed_Range_Global_Mean_Out_Tab),OAC_Converted_Transformed_Range_Global_Mean_Out_Tab)
colnames(OAC_Converted_Transformed_Range_Global_Mean_Out_VarID) <- c("OAC_Converted_Transformed_Range_Global_Mean_VarID",names(OAC_Converted_Transformed_Range_Global_Mean_Out_Tab))

OAC_Converted_Transformed_Range_Global_Mean_Zero_Line<- data.frame(OAC_Converted_Transformed_Range_Global_Mean_Out_VarID$OAC_Converted_Transformed_Range_Global_Mean_VarID, 0)
names(OAC_Converted_Transformed_Range_Global_Mean_Zero_Line)<- c("OAC_Converted_Transformed_Range_Global_Mean_VarID", "OAC_Converted_Transformed_Range_Global_Mean_Zero")

for (nn in 1:CN) {assign(paste("RPLOT_G", nn, sep=""), c(0))}

RPLOTList <-list()

for (ggg in 1:OAC_Converted_Transformed_Range_Global_Mean_Clus_Nu)
{
jjj<-ggg+2
HJ<- paste(BARPLOT_V, ggg, sep=" ")


	if(RQCLUS=="SUPERGROUP"){RADIALTITLE <- paste(HJ, sep="")}
	if(RQCLUS=="GROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, letters[ggg], sep="")
	RADIALTITLE <- paste(BARPLOT_V,BarPlotClus, sep=" ")	
	}
	if(RQCLUS=="SUBGROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, ggg, sep="")
	RADIALTITLE <- paste(BARPLOT_V,BarPlotClus, sep=" ")
	}

RADIALSUB_2 <-paste(RADIALSUB_B," - ",BARSUB_G, sep="")

if(RQCLUS=="GROUP" | RQCLUS=="SUBGROUP")
{
RADIALSUB_2 <-paste(RADIALSUB_B," - ",BARSUB_G, sep="")
}

par(cex.main=1.5)

radial.plot(OAC_Converted_Transformed_Range_Global_Mean_Zero_Line$OAC_Converted_Transformed_Range_Global_Mean_Zero, labels=c(OAC_Converted_Transformed_Range_Global_Mean_Zero_Line$OAC_Converted_Transformed_Range_Global_Mean_VarID), start=0,clockwise=TRUE, rp.type="p", line.col="red",lwd=3, show.grid=TRUE, grid.col="grey", point.symbols=16,point.col="black",show.centroid=FALSE, mar=c(4,4,6,4), radial.lim=range(Plot_Min_Best_G,Plot_Max_Best_G))

radial.plot(OAC_Converted_Transformed_Range_Global_Mean_Out_VarID[,jjj], labels=c(OAC_Converted_Transformed_Range_Global_Mean_Out_VarID$OAC_Converted_Transformed_Range_Global_Mean_VarID), start=0,clockwise=TRUE, rp.type="p&s", line.col=RADIALCOLOUR,lwd=3, show.grid=TRUE, grid.col="grey", point.symbols=16,point.col="black",show.centroid=FALSE, radial.lim=range(Plot_Min_Best_G,Plot_Max_Best_G),  mar=c(4,4,6,4), title(main=get(paste("RADIALTITLE", sep="")), outer=FALSE, line=4), add=T)
mtext(RADIALSUB_2, line=2.5)

assign(paste("RPLOT_G",ggg, sep=""),recordPlot())

RPLOTList <- list(RPLOTList, get(paste("RPLOT_G", ggg, sep="")))

graphics.off()
}

dir.create("Cluster Plots", showWarnings = FALSE)

if(RQOAC =="NO")	
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," Profiles - Radial Plots (Individual Global Mean).pdf", sep=""), title = RADIALNAME_B, family='Courier', width=10, height=10)
}

if(RQOAC =="YES")
{
if(RQCLUS=="SUPERGROUP")
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," Profiles - Radial Plots (Individual Global Mean).pdf", sep=""), title = RADIALNAME_B, family='Courier', width=10, height=10)
}
if(RQCLUS=="GROUP" | RQCLUS=="SUBGROUP")
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," ",Best_Cluster_Input, " Profiles - Radial Plots (Individual Global Mean).pdf", sep=""), title = RADIALNAME_B, family='Courier', width=10, height=10)
}
}
print(RPLOTList)
dev.off()

dir.create("Cluster Plots", showWarnings = FALSE)

if (CN<=8)
{

if(RQOAC =="NO")	
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," Profiles - Radial Plots (Combined Global Mean).pdf", sep=""), title = RADIALNAME_B, family='Courier', width=35, height=20)
}

if(RQOAC =="YES")
{
if(RQCLUS=="SUPERGROUP")
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," Profiles - Radial Plots (Combined Global Mean).pdf", sep=""), title = RADIALNAME_B, family='Courier', width=35, height=20)
}
if(RQCLUS=="GROUP" | RQCLUS=="SUBGROUP")
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," ",Best_Cluster_Input, " Profiles - Radial Plots (Combined Global Mean).pdf", sep=""), title = RADIALNAME_B, family='Courier', width=35, height=20)
}
}

par(mfrow = c(2,4), cex.main=2.5, xpd=TRUE, oma=c(0,0,5,0))

for (ggg in 1:OAC_Converted_Transformed_Range_Global_Mean_Clus_Nu)
{
jjj<-ggg+2

HJ<- paste(BARPLOT_V, ggg, sep=" ")

	if(RQCLUS=="SUPERGROUP"){RADIALTITLE <- paste(HJ, sep="")}
	if(RQCLUS=="GROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, letters[ggg], sep="")
	RADIALTITLE <- paste(BARPLOT_V,BarPlotClus, sep=" ")	
	}
	if(RQCLUS=="SUBGROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, ggg, sep="")
	RADIALTITLE <- paste(BARPLOT_V,BarPlotClus, sep=" ")
	}

RADIALNAME_2 <-paste(RADIALNAME_B," - ",BARSUB_G, sep="")

if(RQCLUS=="GROUP" | RQCLUS=="SUBGROUP")
{
RADIALNAME_2 <-paste(RADIALNAME_B," - ",BARSUB_G, sep="")
}

radial.plot(OAC_Converted_Transformed_Range_Global_Mean_Zero_Line$OAC_Converted_Transformed_Range_Global_Mean_Zero, labels=c(OAC_Converted_Transformed_Range_Global_Mean_Zero_Line$OAC_Converted_Transformed_Range_Global_Mean_VarID), start=0,clockwise=TRUE, rp.type="p", line.col="red",lwd=3, show.grid=TRUE, grid.col="grey", point.symbols=16,point.col="black",show.centroid=FALSE, radial.lim=range(Plot_Min_Best_G,Plot_Max_Best_G), mar=c(5,5,15,5))

radial.plot(OAC_Converted_Transformed_Range_Global_Mean_Out_VarID[,jjj], labels=c(OAC_Converted_Transformed_Range_Global_Mean_Out_VarID$OAC_Converted_Transformed_Range_Global_Mean_VarID), start=0,clockwise=TRUE, rp.type="p&s", line.col=RADIALCOLOUR,lwd=3, show.grid=TRUE, grid.col="grey",
point.symbols=16,point.col="black",show.centroid=FALSE, radial.lim=range(Plot_Min_Best_G,Plot_Max_Best_G), mar=c(5,5,15,5), title(main=get(paste("RADIALTITLE", sep="")), outer=FALSE), add=T)


}
title(main=paste(RADIALNAME_2, sep=""), outer=TRUE, line=-2, cex.main=6)
graphics.off()
}

if (CN>8)
{

if(RQOAC =="NO")	
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," Profiles - Radial Plots (Combined Global Mean).pdf", sep=""), title = RADIALNAME_B, family='Courier', width=35, height=20)
}

if(RQOAC =="YES")
{
if(RQCLUS=="SUPERGROUP")
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," Profiles - Radial Plots (Combined Global Mean).pdf", sep=""), title = RADIALNAME_B, family='Courier', width=35, height=20)
}
if(RQCLUS=="GROUP" | RQCLUS=="SUBGROUP")
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," ",Best_Cluster_Input, " Profiles - Radial Plots (Combined Global Mean).pdf", sep=""), title = RADIALNAME_B, family='Courier', width=35, height=20)
}
}

par(mfrow = c(2,4), cex.main=2.5, xpd=FALSE, mar=c(2,2,2,2))

for (ggg in 1:OAC_Converted_Transformed_Range_Global_Mean_Clus_Nu)
{
jjj<-ggg+2

HJ<- paste(BARPLOT_V, ggg, sep=" ")

	if(RQCLUS=="SUPERGROUP"){RADIALTITLE <- paste(HJ, sep="")}
	if(RQCLUS=="GROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, letters[ggg], sep="")
	RADIALTITLE <- paste(BARPLOT_V,BarPlotClus, sep=" ")	
	}
	if(RQCLUS=="SUBGROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, ggg, sep="")
	RADIALTITLE <- paste(BARPLOT_V,BarPlotClus, sep=" ")
	}
	
radial.plot(OAC_Converted_Transformed_Range_Global_Mean_Zero_Line$OAC_Converted_Transformed_Range_Global_Mean_Zero, labels=c(OAC_Converted_Transformed_Range_Global_Mean_Zero_Line$OAC_Converted_Transformed_Range_Global_Mean_VarID), start=0,clockwise=TRUE, rp.type="p", line.col="red",lwd=3, show.grid=TRUE, grid.col="grey", point.symbols=16,point.col="black",show.centroid=FALSE, radial.lim=range(Plot_Min_Best_G,Plot_Max_Best_G), mar=c(5,5,15,5))

radial.plot(OAC_Converted_Transformed_Range_Global_Mean_Out_VarID[,jjj], labels=c(OAC_Converted_Transformed_Range_Global_Mean_Out_VarID$OAC_Converted_Transformed_Range_Global_Mean_VarID), start=0,clockwise=TRUE, rp.type="p&s", line.col=RADIALCOLOUR,lwd=3, show.grid=TRUE, grid.col="grey",
point.symbols=16,point.col="black",show.centroid=FALSE, radial.lim=range(Plot_Min_Best_G,Plot_Max_Best_G), mar=c(5,5,15,5), title(main=get(paste("RADIALTITLE", sep="")), outer=FALSE), add=T)
}
graphics.off()
}

Plot_Min_Best <- min(Plot_Min_Best_C,Plot_Min_Best_G)
Plot_Max_Best <- max(Plot_Max_Best_C,Plot_Max_Best_G)

for (nn in 1:CN) {assign(paste("RPLOT_C", nn, sep=""), c(0))}

RPLOTList <-list()

for (ggg in 1:OAC_Converted_Transformed_Range_Global_Mean_Clus_Nu)
{
jjj<-ggg+2
HJ<- paste(BARPLOT_V, ggg, sep=" ")


	if(RQCLUS=="SUPERGROUP"){RADIALTITLE <- paste(HJ, sep="")}
	if(RQCLUS=="GROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, letters[ggg], sep="")
	RADIALTITLE <- paste(BARPLOT_V,BarPlotClus, sep=" ")	
	}
	if(RQCLUS=="SUBGROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, ggg, sep="")
	RADIALTITLE <- paste(BARPLOT_V,BarPlotClus, sep=" ")
	}

RADIALSUB_3 <-paste(RADIALSUB_B," - ",BARSUB_CG, sep="")

if(RQCLUS=="GROUP" | RQCLUS=="SUBGROUP")
{
RADIALSUB_3 <-paste(RADIALSUB_B," - ",BARSUB_CG, sep="")
}

par(cex.main=1.5)

radial.plot(OAC_Converted_Transformed_Range_Global_Mean_Zero_Line$OAC_Converted_Transformed_Range_Global_Mean_Zero, labels=c(OAC_Converted_Transformed_Range_Global_Mean_Zero_Line$OAC_Converted_Transformed_Range_Global_Mean_VarID), start=0,clockwise=TRUE, rp.type="p", line.col="red",lwd=3, show.grid=TRUE, grid.col="grey", point.symbols=16,point.col="black",show.centroid=FALSE, mar=c(4,4,6,4), radial.lim=range(Plot_Min_Best,Plot_Max_Best))

radial.plot(OAC_Converted_Transformed_Range_Cluster_Mean_Out_VarID[,jjj], labels=c(OAC_Converted_Transformed_Range_Cluster_Mean_Out_VarID$OAC_Converted_Transformed_Range_Cluster_Mean_VarID), start=0,clockwise=TRUE, rp.type="p&s", line.col=BARCOLOURLOW,lwd=3, show.grid=TRUE, grid.col="grey", point.symbols=16,point.col=BARCOLOURLOW,show.centroid=FALSE, radial.lim=range(Plot_Min_Best,Plot_Max_Best),  mar=c(4,4,6,4), title(main=get(paste("RADIALTITLE", sep="")), outer=FALSE, line=4), add=T)
mtext(RADIALSUB_3, line=2.5)

radial.plot(OAC_Converted_Transformed_Range_Global_Mean_Out_VarID[,jjj], labels=c(OAC_Converted_Transformed_Range_Global_Mean_Out_VarID$OAC_Converted_Transformed_Range_Global_Mean_VarID), start=0,clockwise=TRUE, rp.type="p&s", line.col=RADIALCOLOUR,lwd=3, show.grid=TRUE, grid.col="grey", point.symbols=16,point.col="black",show.centroid=FALSE, radial.lim=range(Plot_Min_Best,Plot_Max_Best),  mar=c(4,4,6,4), title(main=get(paste("RADIALTITLE", sep="")), outer=FALSE, line=4), add=T)
mtext(RADIALSUB_3, line=2.5)

assign(paste("RPLOT_C",ggg, sep=""),recordPlot())

RPLOTList <- list(RPLOTList, get(paste("RPLOT_C", ggg, sep="")))

graphics.off()
}

dir.create("Cluster Plots", showWarnings = FALSE)

if(RQOAC =="NO")	
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," Profiles - Radial Plots (Individual Global and Parent Cluster Mean).pdf", sep=""), title = RADIALNAME_B, family='Courier', width=10, height=10)
}

if(RQOAC =="YES")
{
if(RQCLUS=="SUPERGROUP")
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," Profiles - Radial Plots (Individual Global and Parent Cluster Mean).pdf", sep=""), title = RADIALNAME_B, family='Courier', width=10, height=10)
}
if(RQCLUS=="GROUP" | RQCLUS=="SUBGROUP")
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," ",Best_Cluster_Input, " Profiles - Radial Plots (Individual Global and Parent Cluster Mean).pdf", sep=""), title = RADIALNAME_B, family='Courier', width=10, height=10)
}
}
print(RPLOTList)
graphics.off()

dir.create("Cluster Plots", showWarnings = FALSE)

if (CN<=8)
{

if(RQOAC =="NO")	
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," Profiles - Radial Plots (Combined Global and Parent Cluster Mean).pdf", sep=""), title = RADIALNAME_B, family='Courier', width=35, height=20)
}

if(RQOAC =="YES")
{
if(RQCLUS=="SUPERGROUP")
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," Profiles - Radial Plots (Combined Global and Parent Cluster Mean).pdf", sep=""), title = RADIALNAME_B, family='Courier', width=35, height=20)
}
if(RQCLUS=="GROUP" | RQCLUS=="SUBGROUP")
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," ",Best_Cluster_Input, " Profiles - Radial Plots (Combined Global and Parent Cluster Mean).pdf", sep=""), title = RADIALNAME_B, family='Courier', width=35, height=20)
}
}

par(mfrow = c(2,4), cex.main=2.5, xpd=TRUE, oma=c(0,0,5,0))

for (ggg in 1:OAC_Converted_Transformed_Range_Global_Mean_Clus_Nu)
{
jjj<-ggg+2

HJ<- paste(BARPLOT_V, ggg, sep=" ")

	if(RQCLUS=="SUPERGROUP"){RADIALTITLE <- paste(HJ, sep="")}
	if(RQCLUS=="GROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, letters[ggg], sep="")
	RADIALTITLE <- paste(BARPLOT_V,BarPlotClus, sep=" ")	
	}
	if(RQCLUS=="SUBGROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, ggg, sep="")
	RADIALTITLE <- paste(BARPLOT_V,BarPlotClus, sep=" ")
	}

RADIALNAME_3 <-paste(RADIALNAME_B," - ",BARSUB_CG, sep="")

if(RQCLUS=="GROUP" | RQCLUS=="SUBGROUP")
{
RADIALNAME_3 <-paste(RADIALNAME_B," - ",BARSUB_CG, sep="")
}

radial.plot(OAC_Converted_Transformed_Range_Global_Mean_Zero_Line$OAC_Converted_Transformed_Range_Global_Mean_Zero, labels=c(OAC_Converted_Transformed_Range_Global_Mean_Zero_Line$OAC_Converted_Transformed_Range_Global_Mean_VarID), start=0,clockwise=TRUE, rp.type="p", line.col="red",lwd=3, show.grid=TRUE, grid.col="grey", point.symbols=16,point.col="black",show.centroid=FALSE, radial.lim=range(Plot_Min_Best,Plot_Max_Best), mar=c(5,5,15,5))

radial.plot(OAC_Converted_Transformed_Range_Cluster_Mean_Out_VarID[,jjj], labels=c(OAC_Converted_Transformed_Range_Cluster_Mean_Out_VarID$OAC_Converted_Transformed_Range_Cluster_Mean_VarID), start=0,clockwise=TRUE, rp.type="p&s", line.col=BARCOLOURLOW,lwd=3, show.grid=TRUE, grid.col="grey",
point.symbols=16,point.col=BARCOLOURLOW,show.centroid=FALSE, radial.lim=range(Plot_Min_Best,Plot_Max_Best), mar=c(5,5,15,5), title(main=get(paste("RADIALTITLE", sep="")), outer=FALSE), add=T)

radial.plot(OAC_Converted_Transformed_Range_Global_Mean_Out_VarID[,jjj], labels=c(OAC_Converted_Transformed_Range_Global_Mean_Out_VarID$OAC_Converted_Transformed_Range_Global_Mean_VarID), start=0,clockwise=TRUE, rp.type="p&s", line.col=RADIALCOLOUR,lwd=3, show.grid=TRUE, grid.col="grey",
point.symbols=16,point.col="black",show.centroid=FALSE, radial.lim=range(Plot_Min_Best,Plot_Max_Best), mar=c(5,5,15,5), title(main=get(paste("RADIALTITLE", sep="")), outer=FALSE), add=T)
}
title(main=paste(RADIALNAME_3, sep=""), outer=TRUE, line=-2, cex.main=6)
graphics.off()
}

if (CN>8)
{
if(RQOAC =="NO")	
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," Profiles - Radial Plots (Combined Global and Parent Cluster Mean).pdf", sep=""), title = RADIALNAME_B, family='Courier', width=35, height=20)
}

if(RQOAC =="YES")
{
if(RQCLUS=="SUPERGROUP")
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," Profiles - Radial Plots (Combined Global and Parent Cluster Mean).pdf", sep=""), title = RADIALNAME_B, family='Courier', width=35, height=20)
}
if(RQCLUS=="GROUP" | RQCLUS=="SUBGROUP")
{
pdf(file = paste("Cluster Plots/",DOCFILE_B," ",Best_Cluster_Input, " Profiles - Radial Plots (Combined Global and Parent Cluster Mean).pdf", sep=""), title = RADIALNAME_B, family='Courier', width=35, height=20)
}
}

par(mfrow = c(2,4), cex.main=2.5, xpd=FALSE, mar=c(2,2,2,2))

for (ggg in 1:OAC_Converted_Transformed_Range_Global_Mean_Clus_Nu)
{
jjj<-ggg+2

HJ<- paste(BARPLOT_V, ggg, sep=" ")

	if(RQCLUS=="SUPERGROUP"){RADIALTITLE <- paste(HJ, sep="")}
	if(RQCLUS=="GROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, letters[ggg], sep="")
	RADIALTITLE <- paste(BARPLOT_V,BarPlotClus, sep=" ")	
	}
	if(RQCLUS=="SUBGROUP")
	{
	BarPlotClus<-paste(Best_Cluster_Input, ggg, sep="")
	RADIALTITLE <- paste(BARPLOT_V,BarPlotClus, sep=" ")
	}
	
radial.plot(OAC_Converted_Transformed_Range_Global_Mean_Zero_Line$OAC_Converted_Transformed_Range_Global_Mean_Zero, labels=c(OAC_Converted_Transformed_Range_Global_Mean_Zero_Line$OAC_Converted_Transformed_Range_Global_Mean_VarID), start=0,clockwise=TRUE, rp.type="p", line.col="red",lwd=3, show.grid=TRUE, grid.col="grey", point.symbols=16,point.col="black",show.centroid=FALSE, radial.lim=range(Plot_Min_Best,Plot_Max_Best), mar=c(5,5,15,5))

radial.plot(OAC_Converted_Transformed_Range_Cluster_Mean_Out_VarID[,jjj], labels=c(OAC_Converted_Transformed_Range_Cluster_Mean_Out_VarID$OAC_Converted_Transformed_Range_Cluster_Mean_VarID), start=0,clockwise=TRUE, rp.type="p&s", line.col=BARCOLOURLOW,lwd=3, show.grid=TRUE, grid.col="grey",
point.symbols=16,point.col=BARCOLOURLOW,show.centroid=FALSE, radial.lim=range(Plot_Min_Best,Plot_Max_Best), mar=c(5,5,15,5), title(main=get(paste("RADIALTITLE", sep="")), outer=FALSE), add=T)

radial.plot(OAC_Converted_Transformed_Range_Global_Mean_Out_VarID[,jjj], labels=c(OAC_Converted_Transformed_Range_Global_Mean_Out_VarID$OAC_Converted_Transformed_Range_Global_Mean_VarID), start=0,clockwise=TRUE, rp.type="p&s", line.col=RADIALCOLOUR,lwd=3, show.grid=TRUE, grid.col="grey",
point.symbols=16,point.col="black",show.centroid=FALSE, radial.lim=range(Plot_Min_Best,Plot_Max_Best), mar=c(5,5,15,5), title(main=get(paste("RADIALTITLE", sep="")), outer=FALSE), add=T)
}
graphics.off()
}
}

####################################################################################################
# Produce Outputs per Cluster ######################################################################
####################################################################################################

if(RQCLUS=="SUPERGROUP")
{
for (ttt in 1:CN)
{
assign(paste("OAC_Converted_Transformed_Range_Cluster_",ttt, sep=""), OAC_Converted_Transformed_Range_CSV_Output[which(OAC_Converted_Transformed_Range_CSV_Output$Cluster==ttt),1:c(ncol(OAC_Converted_Transformed_Range_CSV_Output))])
}
}

if(RQCLUS=="GROUP")
{
for (ttt in 1:CN)
{
ClusterNumLet<-paste(Best_Cluster_Input, letters[ttt], sep="")
assign(paste("OAC_Converted_Transformed_Range_Cluster_", ClusterNumLet, sep=""), OAC_Converted_Transformed_Range_CSV_Output[which(OAC_Converted_Transformed_Range_CSV_Output$Cluster==ClusterNumLet),1:c(ncol(OAC_Converted_Transformed_Range_CSV_Output))])
}
}

if(RQCLUS=="SUBGROUP")
{
for (ttt in 1:CN)
{
ClusterNumLet<-paste(Best_Cluster_Input, ttt, sep="")
assign(paste("OAC_Converted_Transformed_Range_Cluster_", ClusterNumLet, sep=""), OAC_Converted_Transformed_Range_CSV_Output[which(OAC_Converted_Transformed_Range_CSV_Output$Cluster==ClusterNumLet),1:c(ncol(OAC_Converted_Transformed_Range_CSV_Output))])
}
}


####################################################################################################
# Variable Numbers and Names #######################################################################
####################################################################################################

OAC_Converted_Transformed_Range_Variable_Codes<- data.frame(names(OAC_Converted_Transformed_Range[,1:c(ncol(OAC_Converted_Transformed_Range)-1)]))
OAC_Converted_Transformed_Range_Variable_Codes_and_Numbers<-data.frame (1:nrow(OAC_Converted_Transformed_Range_Variable_Codes),OAC_Converted_Transformed_Range_Variable_Codes)
names(OAC_Converted_Transformed_Range_Variable_Codes_and_Numbers)<- c("Variable Number", "Variable Code")

if(exists("OAC_Input_Lookup")=="TRUE")
{
OAC_Converted_Transformed_Range_Variable_Name_Lookup<-data.frame(OAC_Input_Lookup$VariableCode,OAC_Input_Lookup$VariableClustered)
OAC_Converted_Transformed_Range_Variable_Codes_and_Numbers <- merge(OAC_Converted_Transformed_Range_Variable_Codes_and_Numbers, OAC_Converted_Transformed_Range_Variable_Name_Lookup, by.x = "Variable Code", by.y = "OAC_Input_Lookup.VariableCode")
OAC_Converted_Transformed_Range_Variable_Codes_and_Numbers <- OAC_Converted_Transformed_Range_Variable_Codes_and_Numbers[,c(2,1,3)]
colnames(OAC_Converted_Transformed_Range_Variable_Codes_and_Numbers) <-c("Variable Number", "Variable Code","Variable Description")
}

####################################################################################################
# CSV Outputs ######################################################################################
####################################################################################################

#VARIABLE NUMBERS AND NAMES

dir.create("Variable Metadata", showWarnings = FALSE)

write.table(OAC_Converted_Transformed_Range_Variable_Codes_and_Numbers, paste("Variable Metadata/OAC_Converted_Transformed_Range_Variable_Lookup.csv"), sep = ",", row.names=FALSE, col.names = TRUE, qmethod = "double")

#INDIVIDUAL CLUSTERS

dir.create("Cluster Data/Individual Clusters", showWarnings = FALSE)

if(RQCLUS=="SUPERGROUP")
{
for (ttt in 1:CN)
{
write.table(get(paste("OAC_Converted_Transformed_Range_Cluster_",ttt, sep="")), paste("Cluster Data/Individual Clusters/OAC_Converted_Transformed_Range_", KM, "_KMeans_Runs_Cluster_",ttt,".csv", sep = ""), sep = ",", row.names=FALSE, col.names = TRUE, qmethod = "double")
}
}

if(RQCLUS=="GROUP")
{
for (ttt in 1:CN)
{
ClusterNumLet<-paste(Best_Cluster_Input, letters[ttt], sep="")
write.table(get(paste("OAC_Converted_Transformed_Range_Cluster_",ClusterNumLet, sep="")), paste("Cluster Data/Individual Clusters/OAC_Converted_Transformed_Range_", KM, "_KMeans_Runs_Cluster_",ClusterNumLet,".csv", sep = ""), sep = ",", row.names=FALSE, col.names = TRUE, qmethod = "double")
}
}

if(RQCLUS=="SUBGROUP")
{
for (ttt in 1:CN)
{
ClusterNumLet<-paste(Best_Cluster_Input, ttt, sep="")
write.table(get(paste("OAC_Converted_Transformed_Range_Cluster_",ClusterNumLet, sep="")), paste("Cluster Data/Individual Clusters/OAC_Converted_Transformed_Range_", KM, "_KMeans_Runs_Cluster_",ClusterNumLet,".csv", sep = ""), sep = ",", row.names=FALSE, col.names = TRUE, qmethod = "double")
}
}

#CLUSTERING STATISTICS

dir.create("Cluster Data", showWarnings = FALSE)

write.table(OAC_Converted_Transformed_Range_Best_WSS_A, paste("Cluster Data/OAC_Converted_Transformed_Range_Sum_of_Squares_Metadata_", KM, "_Runs.csv", sep = ""), sep = ",", row.names= FALSE, col.names = TRUE, qmethod = "double")

####################################################################################################
# Cluster SED Calculations #########################################################################
####################################################################################################

CNX <-CN+1
OAC_SED_Metadata <- OAC_Converted_Transformed_Range_Cluster_Metadata[,-(1:4),drop=FALSE] 
OAC_SED_Metadata <- OAC_SED_Metadata[-CNX,,drop=FALSE]

EUC_DIST = function(x1, x2)
{
   temp = x1-x2
   sum(temp*temp)
}

if (RQSED =="YES")
{
for (n in 1:CN) {assign(paste("DistArray", n, sep=""), c(0))}

DistArraySum<-c(0)
DistArrayFinal<-c()
DistArray <- c()
 
pb2 <- txtProgressBar(min = 0, max = RowNumber, style = 3)
 
for(j in 1:RowNumber)
{
  for(p in 1:ColNumber)
  {
  	assign(paste("b", p, sep=""),c(0))
  }
  		
 for(k in 1:CN)
 {	
	for(p in 1:ColNumber)
  	{
  	assign(paste("b", p, sep=""),EUC_DIST(OAC_Converted_Transformed_Range[,p][j],OAC_SED_Metadata[,p][k]))
	}
	Best_Sum<-0
	for (c in 1:ColNumber)
	{
  	Best_Sum <- sum(get(paste("b",c, sep="")),Best_Sum)
  	}
	Best_Sum <-sqrt(Best_Sum)
 	DistArraySum[j] <-Best_Sum
   	DistArrayFinal <-rbind(DistArrayFinal, DistArraySum[j])      
	}
	Sys.sleep(0.1)
	setTxtProgressBar(pb2, j)
	#Progress Bar will appear below. (n.b. It may take some time to appear for larger calculations)
}

   DistArray<-matrix(DistArrayFinal, ncol=k, byrow=TRUE)
   DistArray<-data.frame(DistArray)
   
for (ee in 1:CN)
		{
		if(RQOAC=="NO")	
		{
		if(RQCLUS=="SUPERGROUP")
		{
		names(DistArray)[ee]<-paste("Cluster_", ee,"_SED", sep="")
		}
		if(RQCLUS=="GROUP")
		{
		DistArrayNames<-paste(Best_Cluster_Input, letters[ee], sep="")
		names(DistArray)[ee]<-paste("Cluster_", DistArrayNames,"_SED", sep="")
		}
		if(RQCLUS=="SUBGROUP")
		{
		DistArrayNames<-paste(Best_Cluster_Input, ee, sep="")
		names(DistArray)[ee]<-paste("Cluster_", DistArrayNames,"_SED", sep="")
		}
		}
		if(RQOAC=="YES")	
		{
		if(RQCLUS=="SUPERGROUP")
		{
		names(DistArray)[ee]<-paste("Supergroup_", ee,"_SED", sep="")
		}
		if(RQCLUS=="GROUP")
		{
		DistArrayNames<-paste(Best_Cluster_Input, letters[ee], sep="")
		names(DistArray)[ee]<-paste("Group_", DistArrayNames,"_SED", sep="")
		}
		if(RQCLUS=="SUBGROUP")
		{
		DistArrayNames<-paste(Best_Cluster_Input, ee, sep="")
		names(DistArray)[ee]<-paste("Subgroup_", DistArrayNames,"_SED", sep="")
		}
		}
		}

OAC_Converted_Transformed_Range_DistArray_Min <-as.data.frame(apply(DistArray, 1, min))
names(OAC_Converted_Transformed_Range_DistArray_Min)<- c("Min_OA_Best_SED_Value")

# Saving R Data 
save.image("OAC_Converted_Transformed_Range_Clustered_with_SED_Values.RData")
}

####################################################################################################
# Cluster SED Tables ###############################################################################
####################################################################################################

if (RQSED =="YES")
{
OAC_Converted_Transformed_Range_SED_Best_Cluster_Data<-c(as.data.frame(OA),as.data.frame(OAC_Converted_Transformed_Range_DistArray_Min), as.data.frame(Cluster))

names(OAC_Converted_Transformed_Range_SED_Best_Cluster_Data)<- c("OA", "Min_SED_Value", "Cluster")

OAC_Converted_Transformed_Range_SED_Best_Cluster_Data<-c(as.data.frame(OAC_Converted_Transformed_Range_SED_Best_Cluster_Data), as.data.frame(DistArray))

OAC_Converted_Transformed_Range_SED_Best_Cluster_Data <-data.frame(OAC_Converted_Transformed_Range_SED_Best_Cluster_Data)

OAC_Converted_Transformed_Range_SED_Best_Cluster_Data_No_Min<-c(as.data.frame(OA), as.data.frame(Cluster))

names(OAC_Converted_Transformed_Range_SED_Best_Cluster_Data_No_Min)<- c("OA", "Cluster")

OAC_Converted_Transformed_Range_SED_Best_Cluster_Data_No_Min<-c(as.data.frame(OAC_Converted_Transformed_Range_SED_Best_Cluster_Data_No_Min), as.data.frame(DistArray))

OAC_Converted_Transformed_Range_SED_Best_Cluster_Data_No_Min <-data.frame(OAC_Converted_Transformed_Range_SED_Best_Cluster_Data_No_Min)
}

####################################################################################################
# SED Cluster Sums #################################################################################
####################################################################################################

if(RQSED =="YES")
{

OAC_Converted_Transformed_Range_Best_SED_Order<-data.frame(OAC_Converted_Transformed_Range_SED_Best_Cluster_Data_No_Min)

for (vv in 1:CN) {assign(paste("s_B", vv, sep=""), sum(OAC_Converted_Transformed_Range_Best_SED_Order[vv+2]))}

OAC_Converted_Transformed_Range_Best_SED_Sum <-c(0)
OAC_Converted_Transformed_Range_Best_SED_Cluster <-c(0)
OAC_Converted_Transformed_Range_Best_SED_Sum_Total <-c()

for (dd in 1:CN)
{
	OAC_Converted_Transformed_Range_Best_SED_Sum[dd] <- get(paste("s_B", dd, sep=""))
	OAC_Converted_Transformed_Range_Best_SED_Cluster[dd] <-dd
	if(RQCLUS=="GROUP")
	{
	OAC_Converted_Transformed_Range_Best_SED_Cluster[dd]<-paste(Best_Cluster_Input, letters[dd], sep="")
	}
	if(RQCLUS=="SUBGROUP")
	{
	OAC_Converted_Transformed_Range_Best_SED_Cluster[dd]<-paste(Best_Cluster_Input, dd, sep="")
	}
}

OAC_Converted_Transformed_Range_Best_SED_Sum_Total <- cbind(OAC_Converted_Transformed_Range_Best_SED_Cluster, OAC_Converted_Transformed_Range_Best_SED_Sum)

OAC_Converted_Transformed_Range_Best_SED_Sum_Total<-data.frame(OAC_Converted_Transformed_Range_Best_SED_Sum_Total)

names(OAC_Converted_Transformed_Range_Best_SED_Sum_Total)<- c("Cluster_Number", "SED Value")

OAC_Converted_Transformed_Range_Cluster_SED_Metadata <-data.frame(Within_Cluster_Sum_of_Squares, Points_Within_Cluster, BestWSS$centers)
OAC_Converted_Transformed_Range_Cluster_SED_Metadata$Mean_Cluster_Within_Sum_of_Squares<-OAC_Converted_Transformed_Range_Cluster_SED_Metadata[,1]/OAC_Converted_Transformed_Range_Cluster_SED_Metadata[,2]
OAC_Converted_Transformed_Range_Cluster_SED_Metadata<-rbind(OAC_Converted_Transformed_Range_Cluster_SED_Metadata, sapply(OAC_Converted_Transformed_Range_Cluster_SED_Metadata, mean))
OAC_Converted_Transformed_Range_Cluster_SED_Metadata_Cluster<-mixedsort(unique(OAC_Converted_Transformed_Range$Cluster))
OAC_Converted_Transformed_Range_Cluster_SED_Metadata_Cluster<-as.matrix(OAC_Converted_Transformed_Range_Cluster_SED_Metadata_Cluster)
OAC_Converted_Transformed_Range_Cluster_SED_Metadata$Cluster<-rbind(OAC_Converted_Transformed_Range_Cluster_SED_Metadata_Cluster,"Mean")
OAC_Converted_Transformed_Range_Cluster_SED_MetadataCol<-ncol(OAC_Converted_Transformed_Range_Cluster_SED_Metadata)
OAC_Converted_Transformed_Range_Cluster_SED_MetadataMean<-OAC_Converted_Transformed_Range_Cluster_SED_MetadataCol-1
OAC_Converted_Transformed_Range_Clusters_LastVar<-OAC_Converted_Transformed_Range_Cluster_SED_MetadataCol-2
OAC_Converted_Transformed_Range_Cluster_SED_Metadata<-OAC_Converted_Transformed_Range_Cluster_SED_Metadata[,c(OAC_Converted_Transformed_Range_Cluster_SED_MetadataCol,1:2, OAC_Converted_Transformed_Range_Cluster_SED_MetadataMean, 3:OAC_Converted_Transformed_Range_Clusters_LastVar)] 
OAC_Converted_Transformed_Range_Cluster_SED_Metadata[1]<-c(OAC_Converted_Transformed_Range_Best_SED_Cluster,"Mean")
}

####################################################################################################
# Cluster SED Sums and Averages ####################################################################
####################################################################################################

if (RQSED == "YES")
{
SED_Subset_Combined <-c()
OAC_Converted_Transformed_Range_SED_Subset_Sum_and_Averages <-c()

for (ccc in 1:CN)
	{
	OAC_Converted_Transformed_Range_SED_Subset<- OAC_Converted_Transformed_Range_SED_Best_Cluster_Data[2][OAC_Converted_Transformed_Range_SED_Best_Cluster_Data[3] == ccc]
	SED_Subset_Cluster <-ccc
	if(RQCLUS=="GROUP")
	{
	ClusterSortOrdered<-paste(Best_Cluster_Input, letters[ccc], sep="")
	OAC_Converted_Transformed_Range_SED_Subset<- OAC_Converted_Transformed_Range_SED_Best_Cluster_Data[2][OAC_Converted_Transformed_Range_SED_Best_Cluster_Data[3] == ClusterSortOrdered]
	SED_Subset_Cluster <-ClusterSortOrdered
	}
	if(RQCLUS=="SUBGROUP")
	{
	ClusterSortOrdered<-paste(Best_Cluster_Input, ccc, sep="")
	OAC_Converted_Transformed_Range_SED_Subset<- OAC_Converted_Transformed_Range_SED_Best_Cluster_Data[2][OAC_Converted_Transformed_Range_SED_Best_Cluster_Data[3] == ClusterSortOrdered]
	SED_Subset_Cluster <-ClusterSortOrdered
	}
	SED_Subset_Sum<-sum(OAC_Converted_Transformed_Range_SED_Subset)
	SED_Subset_Mean<-mean(OAC_Converted_Transformed_Range_SED_Subset)
	SED_Subset_Median<-median(OAC_Converted_Transformed_Range_SED_Subset)
	SED_Subset_Mode<-mfv(OAC_Converted_Transformed_Range_SED_Subset)[1]
	SED_Subset_SD<-sd(OAC_Converted_Transformed_Range_SED_Subset)
	SED_Subset_Values <-length(OAC_Converted_Transformed_Range_SED_Subset)
	SED_Subset_Percentages <-(SED_Subset_Values/length(OAC_Converted_Transformed_Range_SED_Best_Cluster_Data[,1]))*100
	SED_Subset_Combined <-cbind(SED_Subset_Cluster, SED_Subset_Values, SED_Subset_Percentages, SED_Subset_Sum, SED_Subset_Mean, SED_Subset_Median, SED_Subset_Mode, SED_Subset_SD)
	OAC_Converted_Transformed_Range_SED_Subset_Sum_and_Averages <- rbind(OAC_Converted_Transformed_Range_SED_Subset_Sum_and_Averages, SED_Subset_Combined)
	}

OAC_Converted_Transformed_Range_SED_Subset_Sum_and_Averages <-data.frame(OAC_Converted_Transformed_Range_SED_Subset_Sum_and_Averages)
names(OAC_Converted_Transformed_Range_SED_Subset_Sum_and_Averages) <- c("Cluster", "Number_Assigned", "Percentage_Assigned", "Cluster_Sum", "Cluster_Mean", "Cluster_Median", "Cluster_Mode", "Cluster_SD")

SED_Subset_Sum_and_Averages<-as.matrix(OAC_Converted_Transformed_Range_SED_Subset_Sum_and_Averages[2:ncol(OAC_Converted_Transformed_Range_SED_Subset_Sum_and_Averages)])
SED_Subset_Sum_and_Averages<-matrix(as.numeric(SED_Subset_Sum_and_Averages),ncol=ncol(SED_Subset_Sum_and_Averages),nrow=nrow(SED_Subset_Sum_and_Averages))
SED_Subset_Sum_and_Averages<-as.data.frame(SED_Subset_Sum_and_Averages)
SED_Subset_Sum_and_AveragesMean<-sapply(SED_Subset_Sum_and_Averages, mean)
SED_Subset_Sum_and_AveragesMeanName<-"Mean"
SED_Subset_Sum_and_AveragesFinal<-as.data.frame(t(c(SED_Subset_Sum_and_AveragesMeanName,SED_Subset_Sum_and_AveragesMean)))

colnames(SED_Subset_Sum_and_AveragesFinal)<-colnames(OAC_Converted_Transformed_Range_SED_Subset_Sum_and_Averages)

OAC_Converted_Transformed_Range_SED_Subset_Sum_and_Averages<-rbind(OAC_Converted_Transformed_Range_SED_Subset_Sum_and_Averages,SED_Subset_Sum_and_AveragesFinal)

SED_Combined <-c()
OAC_Converted_Transformed_Range_SED_Sum_and_Averages <-c()
OAC_Converted_Transformed_Range_SED_Best_Cluster_Data_Only <-(OAC_Converted_Transformed_Range_SED_Best_Cluster_Data[,4:ncol(OAC_Converted_Transformed_Range_SED_Best_Cluster_Data)])

pb8 <- txtProgressBar(min = 0, max = CN, style = 3) 

for (bbb in 1:CN)
	{
	SED_Col <-OAC_Converted_Transformed_Range_SED_Best_Cluster_Data_Only[bbb]
	SED_Col<- apply(SED_Col, 1,as.numeric)
	SED_Sum<-sum(SED_Col)
	SED_Mean<-mean(SED_Col)
	SED_Median<-median(SED_Col)
	SED_Mode<-mfv(SED_Col)[1]
	SED_SD<-sd(SED_Col)
	SED_Cluster <-bbb
	if(RQCLUS=="GROUP")
	{
	SED_Cluster <-paste(Best_Cluster_Input, letters[bbb], sep="")
	}
	if(RQCLUS=="SUBGROUP")
	{
	SED_Cluster <-paste(Best_Cluster_Input, bbb, sep="")
	}
	SED_Values <-length(SED_Col)
	SED_Combined <-cbind(SED_Cluster, SED_Values, SED_Sum, SED_Mean, SED_Median, SED_Mode, SED_SD)
	OAC_Converted_Transformed_Range_SED_Sum_and_Averages <- rbind(OAC_Converted_Transformed_Range_SED_Sum_and_Averages, SED_Combined)
	setTxtProgressBar(pb8, bbb)
	}

OAC_Converted_Transformed_Range_SED_Sum_and_Averages <-data.frame(OAC_Converted_Transformed_Range_SED_Sum_and_Averages)
names(OAC_Converted_Transformed_Range_SED_Sum_and_Averages) <- c("Cluster", "Number_Assigned", "Total_Cluster_Sum", "Total_Cluster_Mean", "Total_Cluster_Median", "Total_Cluster_Mode", "Total_Cluster_SD")

SED_Sum_and_Averages<-as.matrix(OAC_Converted_Transformed_Range_SED_Sum_and_Averages[2:ncol(OAC_Converted_Transformed_Range_SED_Sum_and_Averages)])
SED_Sum_and_Averages<-matrix(as.numeric(SED_Sum_and_Averages),ncol=ncol(SED_Sum_and_Averages),nrow=nrow(SED_Sum_and_Averages))
SED_Sum_and_Averages<-as.data.frame(SED_Sum_and_Averages)
SED_Sum_and_AveragesMean<-sapply(SED_Sum_and_Averages, mean)
SED_Sum_and_AveragesMeanName<-"Mean"
SED_Sum_and_AveragesFinal<-as.data.frame(t(c(SED_Sum_and_AveragesMeanName,SED_Sum_and_AveragesMean)))

colnames(SED_Sum_and_AveragesFinal)<-colnames(OAC_Converted_Transformed_Range_SED_Sum_and_Averages)

OAC_Converted_Transformed_Range_SED_Sum_and_Averages<-rbind(OAC_Converted_Transformed_Range_SED_Sum_and_Averages,SED_Sum_and_AveragesFinal)

}

####################################################################################################
# SED CSV Outputs ##################################################################################
####################################################################################################

if(RQSED =="YES")
{
dir.create("Cluster Data", showWarnings = FALSE)

write.table(OAC_Converted_Transformed_Range_SED_Sum_and_Averages, paste("Cluster Data/OAC_Converted_Transformed_Range_SED_Total_Sum_and_Averages_for_", CN, "_Clusters.csv", sep = ""), sep = ",", row.names= FALSE, col.names = TRUE, qmethod = "double")

write.table(OAC_Converted_Transformed_Range_SED_Subset_Sum_and_Averages, paste("Cluster Data/OAC_Converted_Transformed_Range_SED_Cluster_Sum_and_Averages_for_", CN, "_Clusters.csv", sep = ""), sep = ",", row.names= FALSE, col.names = TRUE, qmethod = "double")

write.csv(OAC_Converted_Transformed_Range_SED_Best_Cluster_Data,file="Cluster Data/OAC_Converted_Transformed_Range_SED_Cluster_Data.csv", row.names= FALSE)
}

####################################################################################################
# Discriminant Projection Cluster Plots ############################################################
####################################################################################################

dir.create("Cluster Plots", showWarnings = FALSE)

DisPlotCol <- (round_any(CN,2, f=ceiling))/2

pdf(file = paste("Cluster Plots/Discriminant Projection Cluster Plots.pdf", sep=""), title = "Discriminant Projection Cluster Plots", family='Courier', width=20, height=20)

palette(rainbow(CN, alpha=0.5))
layout(rbind(1,2), heights=c(7.0,1.0))
par(mar=c(5, 10, 5, 10))
plotcluster(OAC_Converted_Transformed_Range_Input, ClusterPlots, bw= FALSE, pointsbyclvecd = TRUE, pch=16, col= ClusterPlots, main="Discriminant Projection Cluster Plot", xlab = "Discriminant Coordinates", ylab = "Discriminant Coordinates", cex=2.5, cex.main=3, cex.lab = 2.5, cex.axis = 2.5, mgp = c(4,1.5,0))
par(mar=c(0, 0, 0, 0))
plot.new()
if(RQOAC=="NO")
{
if(RQCLUS=="SUPERGROUP")
{	
legend(x="center",y="center", legend = paste("Cluster", 1:CN), title=expression(bold("Clusters")), pch=16, col=1:CN, ncol=DisPlotCol ,bty ="n", cex=2.5, pt.cex=2.5)
}
if(RQCLUS=="GROUP")
{	
legend(x="center",y="center", legend = paste("Cluster ", Best_Cluster_Input, letters[1:CN], sep=""), title=expression(bold("Clusters")), pch=16, col=1:CN, ncol=DisPlotCol ,bty ="n", cex=2.5, pt.cex=2.5)
}
if(RQCLUS=="SUBGROUP")
{	
legend(x="center",y="center", legend = paste("Cluster ", Best_Cluster_Input, 1:CN, sep=""), title=expression(bold("Clusters")), pch=16, col=1:CN, ncol=DisPlotCol ,bty ="n", cex=2.5, pt.cex=2.5)
}
}
if(RQOAC=="YES")
{
if(RQCLUS=="SUPERGROUP")
{	
legend(x="center",y="center", legend = paste("Supergroup", 1:CN), title=expression(bold("Clusters")), pch=16, col=1:CN, ncol=DisPlotCol ,bty ="n", cex=2.5, pt.cex=2.5)
}
if(RQCLUS=="GROUP")
{	
legend(x="center",y="center", legend = paste("Group ", Best_Cluster_Input, letters[1:CN], sep=""), title=expression(bold("Clusters")), pch=16, col=1:CN, ncol=DisPlotCol ,bty ="n", cex=2.5, pt.cex=2.5)
}
if(RQCLUS=="SUBGROUP")
{	
legend(x="center",y="center", legend = paste("Subgroup ", Best_Cluster_Input, 1:CN, sep=""), title=expression(bold("Clusters")), pch=16, col=1:CN, ncol=DisPlotCol ,bty ="n", cex=2.5, pt.cex=2.5)
}
}
assign(paste("DisPlot1", sep=""),recordPlot())

palette(rainbow(CN, alpha=0.5))
layout(rbind(1,2), heights=c(7.0,1.0))
par(mar=c(5, 10, 5, 10))
plotcluster(OAC_Converted_Transformed_Range_Input, ClusterPlots, bw= FALSE, pointsbyclvecd = FALSE, col= ClusterPlots, main="Discriminant Projection Cluster Plot", xlab = "Discriminant Coordinates", ylab = "Discriminant Coordinates", cex=2.5, cex.main=3, cex.lab = 2.5, cex.axis = 2.5, mgp = c(4,1.5,0))
par(mar=c(0, 0, 0, 0))
plot.new()
if(RQOAC=="NO")
{
if(RQCLUS=="SUPERGROUP")
{
legend('center',legend = paste("Cluster", 1:CN), title=expression(bold("Clusters")), pch=1:CN, col=1:CN, ncol=DisPlotCol ,bty ="n", cex=2.5, pt.cex=2.5)
}
if(RQCLUS=="GROUP")
{
legend('center',legend = paste("Cluster ", Best_Cluster_Input, letters[1:CN], sep=""), title=expression(bold("Clusters")), pch=1:CN, col=1:CN, ncol=DisPlotCol ,bty ="n", cex=2.5, pt.cex=2.5)
}
if(RQCLUS=="SUBGROUP")
{
legend('center',legend = paste("Cluster ", Best_Cluster_Input, 1:CN, sep=""), title=expression(bold("Clusters")), pch=1:CN, col=1:CN, ncol=DisPlotCol ,bty ="n", cex=2.5, pt.cex=2.5)
}
}
if(RQOAC=="YES")
{
if(RQCLUS=="SUPERGROUP")
{
legend('center',legend = paste("Supergroup", 1:CN), title=expression(bold("Clusters")), pch=1:CN, col=1:CN, ncol=DisPlotCol ,bty ="n", cex=2.5, pt.cex=2.5)
}
if(RQCLUS=="GROUP")
{
legend('center',legend = paste("Group ", Best_Cluster_Input, letters[1:CN], sep=""), title=expression(bold("Clusters")), pch=1:CN, col=1:CN, ncol=DisPlotCol ,bty ="n", cex=2.5, pt.cex=2.5)
}
if(RQCLUS=="SUBGROUP")
{
legend('center',legend = paste("Subgroup ", Best_Cluster_Input, 1:CN, sep=""), title=expression(bold("Clusters")), pch=1:CN, col=1:CN, ncol=DisPlotCol ,bty ="n", cex=2.5, pt.cex=2.5)
}
}
assign(paste("DisPlot2", sep=""),recordPlot())

palette("default")
blackpalette <- adjustcolor(palette(), alpha.f = 0.5)
palette(blackpalette)
layout(rbind(1,2), heights=c(7.0,1.0))
par(mar=c(5, 10, 5, 10))
plotcluster(OAC_Converted_Transformed_Range_Input, ClusterPlots, bw= TRUE, pointsbyclvecd = FALSE, main="Discriminant Projection Cluster Plot", xlab = "Discriminant Coordinates", ylab = "Discriminant Coordinates", cex=2.5, cex.main=3, cex.lab = 2.5, cex.axis = 2.5, mgp = c(4,1.5,0))
par(mar=c(0, 0, 0, 0))
plot.new()
if(RQOAC=="NO")
{
if(RQCLUS=="SUPERGROUP")
{
legend('center',legend = paste("Cluster", 1:CN), title=expression(bold("Clusters")), pch=1:CN, ncol=DisPlotCol ,bty ="n", cex=2.5, pt.cex=2.5)
}
if(RQCLUS=="GROUP")
{
legend('center',legend = paste("Cluster ", Best_Cluster_Input, letters[1:CN], sep=""), title=expression(bold("Clusters")), pch=1:CN, ncol=DisPlotCol ,bty ="n", cex=2.5, pt.cex=2.5)
}
if(RQCLUS=="SUBGROUP")
{
legend('center',legend = paste("Cluster ", Best_Cluster_Input, 1:CN, sep=""), title=expression(bold("Clusters")), pch=1:CN, ncol=DisPlotCol ,bty ="n", cex=2.5, pt.cex=2.5)
}
}
if(RQOAC=="YES")
{
if(RQCLUS=="SUPERGROUP")
{
legend('center',legend = paste("Supergroup", 1:CN), title=expression(bold("Clusters")), pch=1:CN, ncol=DisPlotCol ,bty ="n", cex=2.5, pt.cex=2.5)
}
if(RQCLUS=="GROUP")
{
legend('center',legend = paste("Group ", Best_Cluster_Input, letters[1:CN], sep=""), title=expression(bold("Clusters")), pch=1:CN, ncol=DisPlotCol ,bty ="n", cex=2.5, pt.cex=2.5)
}
if(RQCLUS=="SUBGROUP")
{
legend('center',legend = paste("Subgroup ", Best_Cluster_Input, 1:CN, sep=""), title=expression(bold("Clusters")), pch=1:CN, ncol=DisPlotCol ,bty ="n", cex=2.5, pt.cex=2.5)
}
}
assign(paste("DisPlot3", sep=""),recordPlot())

graphics.off()

if (Sys.info()['sysname']=="Windows")	
{
windowsFonts(Courier=windowsFont("TT Courier New"))
png(filename = paste("Cluster Plots/Discriminant Projection Cluster Plot.png", sep=""), width = 50, height = 50, units = "cm", res=600, family = 'Courier')
} else
{
png(filename = paste("Cluster Plots/Discriminant Projection Cluster Plot.png", sep=""), width = 50, height = 50, units = "cm", res=600, family = 'Courier', type= 'cairo-png')
}
palette(rainbow(CN, alpha=0.5))
layout(rbind(1,2), heights=c(7.0,1.0))
par(mar=c(5, 10, 5, 10))
plotcluster(OAC_Converted_Transformed_Range_Input, ClusterPlots, bw= FALSE, pointsbyclvecd = TRUE, pch=16, col= ClusterPlots, main="Discriminant Projection Cluster Plot", xlab = "Discriminant Coordinates", ylab = "Discriminant Coordinates", cex=2.5, cex.main=3, cex.lab = 2.5, cex.axis = 2.5, mgp = c(4,1.5,0))
par(mar=c(0, 0, 0, 0))
plot.new()
if(RQOAC=="NO")
{
if(RQCLUS=="SUPERGROUP")
{	
legend(x="center",y="center", legend = paste("Cluster", 1:CN), title=expression(bold("Clusters")), pch=16, col=1:CN, ncol=DisPlotCol ,bty ="n", cex=2.5, pt.cex=2.5)
}
if(RQCLUS=="GROUP")
{	
legend(x="center",y="center", legend = paste("Cluster ", Best_Cluster_Input, letters[1:CN], sep=""), title=expression(bold("Clusters")), pch=16, col=1:CN, ncol=DisPlotCol ,bty ="n", cex=2.5, pt.cex=2.5)
}
if(RQCLUS=="SUBGROUP")
{	
legend(x="center",y="center", legend = paste("Cluster ", Best_Cluster_Input, 1:CN, sep=""), title=expression(bold("Clusters")), pch=16, col=1:CN, ncol=DisPlotCol ,bty ="n", cex=2.5, pt.cex=2.5)
}
}
if(RQOAC=="YES")
{
if(RQCLUS=="SUPERGROUP")
{	
legend(x="center",y="center", legend = paste("Supergroup", 1:CN), title=expression(bold("Clusters")), pch=16, col=1:CN, ncol=DisPlotCol ,bty ="n", cex=2.5, pt.cex=2.5)
}
if(RQCLUS=="GROUP")
{	
legend(x="center",y="center", legend = paste("Group ", Best_Cluster_Input, letters[1:CN], sep=""), title=expression(bold("Clusters")), pch=16, col=1:CN, ncol=DisPlotCol ,bty ="n", cex=2.5, pt.cex=2.5)
}
if(RQCLUS=="SUBGROUP")
{	
legend(x="center",y="center", legend = paste("Subgroup ", Best_Cluster_Input, 1:CN, sep=""), title=expression(bold("Clusters")), pch=16, col=1:CN, ncol=DisPlotCol ,bty ="n", cex=2.5, pt.cex=2.5)
}
}

graphics.off()

####################################################################################################
# Lorenz Curves and Gini Coefficients ##############################################################
####################################################################################################

dir.create("Cluster Evaluation", showWarnings = FALSE)

OAC_Converted_Transformed_Range_Matrix<-as.matrix(OAC_Converted_Transformed_Range_Input)

OldRange = max(OAC_Converted_Transformed_Range_Matrix) - min(OAC_Converted_Transformed_Range_Matrix)
NewRange = (1 - 0)
OAC_Converted_Transformed_Range_MatrixS = (((OAC_Converted_Transformed_Range_Matrix - min(OAC_Converted_Transformed_Range_Matrix)) * NewRange) / OldRange) + 0
OAC_Converted_Transformed_Range_Matrix<-OAC_Converted_Transformed_Range_MatrixS

OAC_Converted_Transformed_Range_MatrixVarNum <- ncol(OAC_Converted_Transformed_Range_Matrix)
OAC_Converted_Transformed_Range_LorenzRbind <-list()
OAC_Converted_Transformed_Range_GiniVar <- list()

for (lcc in 1:OAC_Converted_Transformed_Range_MatrixVarNum)
{
vlcc <-paste("s", lcc, sep="")
OAC_Converted_Transformed_Range_MatrixCol<-OAC_Converted_Transformed_Range_Matrix[,lcc]
gini<-ineq(OAC_Converted_Transformed_Range_MatrixCol,type="Gini")
OAC_Converted_Transformed_Range_LorenzCol<-Lc(OAC_Converted_Transformed_Range_MatrixCol)
OAC_Converted_Transformed_Range_LorenzCbind <- cbind(lcc, OAC_Converted_Transformed_Range_LorenzCol$p, OAC_Converted_Transformed_Range_LorenzCol$L)
row.names(OAC_Converted_Transformed_Range_LorenzCbind) <- NULL
OAC_Converted_Transformed_Range_LorenzCbind<-data.frame(OAC_Converted_Transformed_Range_LorenzCbind)
colnames(OAC_Converted_Transformed_Range_LorenzCbind) <-c("Variable", "p", "L")
assign(paste("OAC_Converted_Transformed_Range_LC", lcc, sep=""), OAC_Converted_Transformed_Range_LorenzCbind)
OAC_Converted_Transformed_Range_LorenzRbind[[lcc]] <- data.frame(OAC_Converted_Transformed_Range_LorenzCbind)
GiniVar <- cbind(data.frame(lcc, gini))
colnames(GiniVar)<-c("Variable", "Gini")
OAC_Converted_Transformed_Range_GiniVar[[lcc]] <- data.frame(GiniVar)
}

OAC_Converted_Transformed_Range_LorenzRbindFill <- rbind.fill(OAC_Converted_Transformed_Range_LorenzRbind)
OAC_Converted_Transformed_Range_GiniRbindFill <- rbind.fill(OAC_Converted_Transformed_Range_GiniVar)

NVar <-ncol(OAC_Converted_Transformed_Range_Input)
MinSLC<-0
NVarLists <- round_any((NVar/10), 1, f =floor)

for (slc in 1:NVarLists)
{
MinSLCValue <- MinSLC+1
MaxSLCValue <- MinSLC+10
assign(paste("SLCRange", slc, sep=""), MinSLCValue:MaxSLCValue)
MinSLC<- MaxSLCValue
}
MinSLCValue <- MinSLC+1
MaxSLCValue <- max(NVar)
assign(paste("SLCRange", slc+1, sep=""), MinSLCValue:MaxSLCValue)

NVarAdd=1
if((NVarLists*10)==ColNumber){NVarAdd=0}

NVarListsAll <- NVarLists+NVarAdd
LorenzPlotList <- list()

for (plc in 1:NVarListsAll)
{
plcrange <- get(paste("SLCRange", plc, sep=""))
plcrangeNum <- count(plcrange)
plcrangeNumMax <- sum(plcrangeNum$freq)
OAC_Converted_Transformed_Range_LorenzLorenzPlot <- subset(OAC_Converted_Transformed_Range_LorenzRbindFill, (Variable %in% c(plcrange)))
OAC_Converted_Transformed_Range_LorenzLorenzPlot$Variable<- as.factor(OAC_Converted_Transformed_Range_LorenzLorenzPlot$Variable)
XYLine<-data.frame(p=0:1, L=0:1, Variable="x=y")
OAC_Converted_Transformed_Range_LorenzPlotData<-rbind(XYLine, OAC_Converted_Transformed_Range_LorenzLorenzPlot)
palette(rainbow(plcrangeNumMax, alpha=0.5))

LorenzPlot <- ggplot(data=subset(OAC_Converted_Transformed_Range_LorenzPlotData,Variable!='x=y'), aes(p, L, group = Variable, colour=Variable))
LorenzPlot<- LorenzPlot + geom_path(alpha = 0.8) + ggtitle("Lorenz Curves\n") + theme(plot.title = element_text(size = 40, colour = "black", face = "bold")) + theme(legend.title = element_text(size=34)) + theme(legend.text = element_text(size=30)) + theme(axis.title.x = element_text(size=34,colour = "black"),axis.text.x=element_text(angle=0, vjust=0.0, hjust=0.5, size=30))+theme(axis.title.y = element_text(size=34,colour = "black"),axis.text.y=element_text(size=30)) + guides(colour = guide_legend(title = "Variable", keywidth = 3, keyheight = 3, order = 1, override.aes = list(size = 3.0))) + xlab("\n x") + ylab("y\n") + theme (plot.margin = unit(c(1, 1, 1, 1), "cm")) + theme(legend.key = element_blank())
LorenzPlot<- LorenzPlot + geom_line(data=subset(OAC_Converted_Transformed_Range_LorenzPlotData,Variable =='x=y'), alpha=0.8, aes(p, L, size= ''), colour="black", inherit.aes = FALSE)+ guides (size = guide_legend(title = "x=y", order=2, keywidth = 3, keyheight = 3, override.aes = list(size = 3.0))) + theme(legend.key = element_blank())+ scale_y_continuous(limits=c(0, 1), breaks = seq(0, 1, 0.25))+ scale_x_continuous(limits=c(0, 1), breaks = seq(0, 1, 0.25))
assign(paste("LorenzPlot", plc, sep=""), LorenzPlot)
LorenzPlotList <- list(LorenzPlotList, get(paste("LorenzPlot", plc, sep="")))
}

pdf(file = paste("Cluster Evaluation/Lorenz Curves.pdf", sep=""), title = "Lorenz Curves", family='Courier', width=20, height=20)
LorenzPlotList
graphics.off()

GiniPlot<- ggplot(OAC_Converted_Transformed_Range_GiniRbindFill, aes(x = factor(Variable), y = Gini)) + geom_bar(stat = "identity")+ ggtitle("Gini Coefficients\n") + theme(plot.title = element_text(size = 40, colour = "black", face = "bold")) + theme(axis.title.x = element_text(size=26,colour = "black"),axis.text.x=element_text(angle=270, vjust=0.5, hjust=0.5, size=22)) +theme(axis.title.y = element_text(size=26,colour = "black"),axis.text.y=element_text(size=22)) + xlab("\n Variable Number") + ylab("Gini Coefficient\n")+ scale_y_continuous(limits=c(0, 1), breaks = seq(0, 1, 0.20)) + theme (plot.margin = unit(c(1, 2, 1, 1), "cm"))

pdf(file = paste("Cluster Evaluation/Gini Coefficients.pdf", sep=""), title = paste("Gini Coefficients", sep=""), family='Courier', width=30, height=20)
GiniPlot
graphics.off()

palette("default")

####################################################################################################
# Correlation ######################################################################################
####################################################################################################

OAC_Converted_Transformed_Range_Input_CORR <- data.frame(OAC_Converted_Transformed_Range_Input)

OAC_Converted_Transformed_Range_Input_Variable_Names <-data.frame(colnames(OAC_Converted_Transformed_Range_Input_CORR))
OAC_Converted_Transformed_Range_Input_Variable_Numbers <-nrow(OAC_Converted_Transformed_Range_Input_Variable_Names)
OAC_Converted_Transformed_Range_Input_Variable_Number_Range <-data.frame(1:OAC_Converted_Transformed_Range_Input_Variable_Numbers)
OAC_Converted_Transformed_Range_Input_Variables <- cbind(OAC_Converted_Transformed_Range_Input_Variable_Number_Range, OAC_Converted_Transformed_Range_Input_Variable_Names)
colnames(OAC_Converted_Transformed_Range_Input_Variables) <-c("Variable Number", "Variable Code")
OAC_Converted_Transformed_Range_Input_Variables_Code<-OAC_Converted_Transformed_Range_Input_Variables
if(exists("OAC_Input_Lookup")=="TRUE")
{
OAC_Converted_Transformed_Range_Variable_Name_Lookup<-data.frame(OAC_Input_Lookup$VariableCode,OAC_Input_Lookup$VariableClustered)
OAC_Converted_Transformed_Range_Input_Variables <- merge(OAC_Converted_Transformed_Range_Input_Variables, OAC_Converted_Transformed_Range_Variable_Name_Lookup, by.x = "Variable Code", by.y = "OAC_Input_Lookup.VariableCode")
OAC_Converted_Transformed_Range_Input_Variables <- OAC_Converted_Transformed_Range_Input_Variables[,c(2,1,3)]
colnames(OAC_Converted_Transformed_Range_Input_Variables) <-c("Variable Number", "Variable Code","Variable Description")
}

OAC_Converted_Transformed_Range_P_CORR <- as.matrix(OAC_Converted_Transformed_Range_Input_CORR,rownames.force=TRUE)

OAC_Converted_Transformed_Range_P_CORR <- as.matrix(OAC_Converted_Transformed_Range_P_CORR)

RowNumber <- nrow(OAC_Converted_Transformed_Range_Input_CORR)
ColNumber <- ncol(OAC_Converted_Transformed_Range_Input_CORR)

R <- rcorr(OAC_Converted_Transformed_Range_P_CORR, type="pearson")$r
p <- rcorr(OAC_Converted_Transformed_Range_P_CORR, type="pearson")$P

CORR_STARS <- ifelse(p == 0, "", ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05,
"* ", " "))))

R <- format(round(cbind(rep(-1.11, ncol(OAC_Converted_Transformed_Range_P_CORR)), R), 2))[,-1]

OAC_Converted_Transformed_Range_P_SIG_CORR <- matrix(paste(R, CORR_STARS, sep=""), ncol=ncol(OAC_Converted_Transformed_Range_P_CORR))
diag(OAC_Converted_Transformed_Range_P_SIG_CORR) <- paste(diag(R), " ", sep="")
rownames(OAC_Converted_Transformed_Range_P_SIG_CORR) <- colnames(OAC_Converted_Transformed_Range_P_CORR)
colnames(OAC_Converted_Transformed_Range_P_SIG_CORR) <- paste(colnames(OAC_Converted_Transformed_Range_P_CORR), "", sep="")

OAC_Converted_Transformed_Range_P_CORR_PEARSON_P <- rcorr(OAC_Converted_Transformed_Range_P_CORR, type="pearson")$P 
OAC_Converted_Transformed_Range_P_CORR_PEARSON_P <- format(round(cbind(rep(-1.11, ncol(OAC_Converted_Transformed_Range_P_CORR)), OAC_Converted_Transformed_Range_P_CORR_PEARSON_P), 4))[,-1] 
OAC_Converted_Transformed_Range_P_SIG_PVALUE <- as.data.frame(OAC_Converted_Transformed_Range_P_CORR_PEARSON_P)
OAC_Converted_Transformed_Range_P_SIG_CORR <- as.data.frame(OAC_Converted_Transformed_Range_P_SIG_CORR)

counta=1
countb=2

OAC_Converted_Transformed_Range_P_SEL_PVALUE <-OAC_Converted_Transformed_Range_P_SIG_PVALUE[1, 1:ncol(OAC_Converted_Transformed_Range_P_SIG_PVALUE)]
OAC_Converted_Transformed_Range_P_SEL_FINAL <-OAC_Converted_Transformed_Range_P_SEL_PVALUE
names(OAC_Converted_Transformed_Range_P_SEL_FINAL)<- c(names(OAC_Converted_Transformed_Range_P_SIG_PVALUE))

for(uu in 2:nrow(OAC_Converted_Transformed_Range_P_SIG_PVALUE))
	{
	OAC_Converted_Transformed_Range_P_SEL_PVALUE <-OAC_Converted_Transformed_Range_P_SIG_PVALUE[uu, countb:ncol(OAC_Converted_Transformed_Range_P_SIG_PVALUE)]
	OAC_Converted_Transformed_Range_P_SEL_CORR <-OAC_Converted_Transformed_Range_P_SIG_CORR[uu, 1:counta]

	OAC_Converted_Transformed_Range_P_SEL_COMBINE <- cbind(OAC_Converted_Transformed_Range_P_SEL_CORR, OAC_Converted_Transformed_Range_P_SEL_PVALUE)[1:ColNumber]
	names(OAC_Converted_Transformed_Range_P_SEL_COMBINE)<- c(names(OAC_Converted_Transformed_Range_P_SIG_PVALUE))
	OAC_Converted_Transformed_Range_P_SEL_FINAL <-rbind(OAC_Converted_Transformed_Range_P_SEL_FINAL,OAC_Converted_Transformed_Range_P_SEL_COMBINE)
	names(OAC_Converted_Transformed_Range_P_SEL_FINAL)<- c(names(OAC_Converted_Transformed_Range_P_SIG_PVALUE))
	counta=counta+1
	countb=countb+1
	}

OAC_Converted_Transformed_Range_S_CORR <- as.matrix(OAC_Converted_Transformed_Range_Input_CORR,rownames.force=TRUE)

OAC_Converted_Transformed_Range_S_CORR <- as.matrix(OAC_Converted_Transformed_Range_S_CORR)

#These calculations may take a few minutes to complete
R <- rcorr(OAC_Converted_Transformed_Range_S_CORR, type="spearman")$r
p <- rcorr(OAC_Converted_Transformed_Range_S_CORR, type="spearman")$P

CORR_STARS <- ifelse(p == 0, "", ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05,
"* ", " "))))

R <- format(round(cbind(rep(-1.11, ncol(OAC_Converted_Transformed_Range_S_CORR)), R), 2))[,-1]

OAC_Converted_Transformed_Range_S_SIG_CORR <- matrix(paste(R, CORR_STARS, sep=""), ncol=ncol(OAC_Converted_Transformed_Range_S_CORR))
diag(OAC_Converted_Transformed_Range_S_SIG_CORR) <- paste(diag(R), " ", sep="")
rownames(OAC_Converted_Transformed_Range_S_SIG_CORR) <- colnames(OAC_Converted_Transformed_Range_S_CORR)
colnames(OAC_Converted_Transformed_Range_S_SIG_CORR) <- paste(colnames(OAC_Converted_Transformed_Range_S_CORR), "", sep="")

#These calculations may take a few minutes to complete
OAC_Converted_Transformed_Range_S_CORR_SPEARMAN_P <- rcorr(OAC_Converted_Transformed_Range_S_CORR, type="spearman")$P 
OAC_Converted_Transformed_Range_S_CORR_SPEARMAN_P <- format(round(cbind(rep(-1.11, ncol(OAC_Converted_Transformed_Range_S_CORR)), OAC_Converted_Transformed_Range_S_CORR_SPEARMAN_P), 4))[,-1] 
OAC_Converted_Transformed_Range_S_SIG_PVALUE <- as.data.frame(OAC_Converted_Transformed_Range_S_CORR_SPEARMAN_P)
OAC_Converted_Transformed_Range_S_SIG_CORR <- as.data.frame(OAC_Converted_Transformed_Range_S_SIG_CORR)

counta=1
countb=2

OAC_Converted_Transformed_Range_S_SEL_PVALUE <-OAC_Converted_Transformed_Range_S_SIG_PVALUE[1, 1:ncol(OAC_Converted_Transformed_Range_S_SIG_PVALUE)]
OAC_Converted_Transformed_Range_S_SEL_FINAL <-OAC_Converted_Transformed_Range_S_SEL_PVALUE
names(OAC_Converted_Transformed_Range_S_SEL_FINAL)<- c(names(OAC_Converted_Transformed_Range_S_SIG_PVALUE))

for(uu in 2:nrow(OAC_Converted_Transformed_Range_S_SIG_PVALUE))
	{
	
	OAC_Converted_Transformed_Range_S_SEL_PVALUE <-OAC_Converted_Transformed_Range_S_SIG_PVALUE[uu, countb:ncol(OAC_Converted_Transformed_Range_S_SIG_PVALUE)]
	OAC_Converted_Transformed_Range_S_SEL_CORR <-OAC_Converted_Transformed_Range_S_SIG_CORR[uu, 1:counta]

	OAC_Converted_Transformed_Range_S_SEL_COMBINE <- cbind(OAC_Converted_Transformed_Range_S_SEL_CORR, OAC_Converted_Transformed_Range_S_SEL_PVALUE)[1:ColNumber]
	names(OAC_Converted_Transformed_Range_S_SEL_COMBINE)<- c(names(OAC_Converted_Transformed_Range_S_SIG_PVALUE))
	OAC_Converted_Transformed_Range_S_SEL_FINAL <-rbind(OAC_Converted_Transformed_Range_S_SEL_FINAL,OAC_Converted_Transformed_Range_S_SEL_COMBINE)
	names(OAC_Converted_Transformed_Range_S_SEL_FINAL)<- c(names(OAC_Converted_Transformed_Range_S_SIG_PVALUE))
	counta=counta+1
	countb=countb+1

	}

OAC_Converted_Transformed_Range_CORR <- as.matrix(OAC_Converted_Transformed_Range_Input_CORR,rownames.force=TRUE)

OAC_Converted_Transformed_Range_CORR <- as.matrix(OAC_Converted_Transformed_Range_CORR)

#These calculations may take a few minutes to complete
OAC_Converted_Transformed_Range_CORR_PEARSON <- rcorr(OAC_Converted_Transformed_Range_CORR, type="pearson") 
OAC_Converted_Transformed_Range_CORR_PEARSON_R <- rcorr(OAC_Converted_Transformed_Range_CORR, type="pearson")$r 
OAC_Converted_Transformed_Range_CORR_PEARSON_P <- rcorr(OAC_Converted_Transformed_Range_CORR, type="pearson")$P 
OAC_Converted_Transformed_Range_CORR_PEARSON_R <- format(round(cbind(rep(-1.11, ncol(OAC_Converted_Transformed_Range_CORR)), OAC_Converted_Transformed_Range_CORR_PEARSON_R), 4))[,-1] 
OAC_Converted_Transformed_Range_CORR_PEARSON_P <- format(round(cbind(rep(-1.11, ncol(OAC_Converted_Transformed_Range_CORR)), OAC_Converted_Transformed_Range_CORR_PEARSON_P), 4))[,-1] 

#These calculations may take a few minutes to complete
OAC_Converted_Transformed_Range_CORR_SPEARMAN = rcorr(OAC_Converted_Transformed_Range_CORR, type="spearman")
OAC_Converted_Transformed_Range_CORR_SPEARMAN_R <- rcorr(OAC_Converted_Transformed_Range_CORR, type="spearman")$r 
OAC_Converted_Transformed_Range_CORR_SPEARMAN_P <- rcorr(OAC_Converted_Transformed_Range_CORR, type="spearman")$P 
OAC_Converted_Transformed_Range_CORR_SPEARMAN_R <- format(round(cbind(rep(-1.11, ncol(OAC_Converted_Transformed_Range_CORR)), OAC_Converted_Transformed_Range_CORR_SPEARMAN_R), 4))[,-1] 
OAC_Converted_Transformed_Range_CORR_SPEARMAN_P <- format(round(cbind(rep(-1.11, ncol(OAC_Converted_Transformed_Range_CORR)), OAC_Converted_Transformed_Range_CORR_SPEARMAN_P), 4))[,-1] 

#Plot Type 1
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot <-OAC_Converted_Transformed_Range_P_CORR
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number<-ncol(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot)
colnames(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot) <- c(1:OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor <- cor(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt <- melt(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor)
colnames(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt) <- c("X1", "X2", "Value")
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt$Value<-cut(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt$Value,breaks=c(-1,-0.60,-0.20,0.20,0.6,1),include.lowest=TRUE,label=c("-1.00 to -0.60 ","-0.60 to -0.20 ","-0.20 to 0.20 ","0.20 to 0.60 ","0.60 to 1.00 "))
CorrMatrix.NoPanel <- list(theme(panel.background=element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank())) 
CorrMatrixPlot <- ggplot(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt,aes(x = X1, y = X2))
CorrMatrixPlot <- CorrMatrixPlot + geom_tile(aes(fill=Value), colour="grey15")
CorrMatrixPlot <- CorrMatrixPlot + scale_fill_manual(values = c('-1.00 to -0.60 ' = "#B2182B", '-0.60 to -0.20 ' = "#F09B7A", '-0.20 to 0.20 ' = "#FFFFFF", '0.20 to 0.60 ' = "#87BEDA", '0.60 to 1.00 ' = "#2166AC"), name="Correlation ", drop=FALSE)
CorrMatrixPlot <- CorrMatrixPlot + scale_x_discrete(name="Variable Number", limits=c(1:OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number))
CorrMatrixPlot <- CorrMatrixPlot + scale_y_discrete(name="Variable Number", limits=c(1:OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number)) 
CorrMatrixPlot <- CorrMatrixPlot + theme(axis.title.x = element_text(size=26,colour = "black", face = "bold"),axis.text.x=element_text(angle=0, size=20))
CorrMatrixPlot <- CorrMatrixPlot + theme(axis.title.y = element_text(size=26,colour = "black", face = "bold"),axis.text.y=element_text(size=20))+ CorrMatrix.NoPanel
CorrMatrixPlot <- CorrMatrixPlot + ggtitle("Correlation Matrix ") + theme(plot.title = element_text(size = 40, colour = "black", face = "bold"))
CorrMatrixPlot <- CorrMatrixPlot + coord_equal() + theme(legend.position = "bottom") 
CorrMatrixPlot <- CorrMatrixPlot + theme(legend.title = element_text(size=26)) + theme(legend.text = element_text(size=20)) 
CorrMatrixPlot <- CorrMatrixPlot + theme(axis.title.x=element_text(hjust=0.5, vjust=0.3))
CorrMatrixPlot <- CorrMatrixPlot + theme(axis.title.y=element_text(hjust=0.5, vjust=0.4))
CorrMatrixPlot <- CorrMatrixPlot + theme(plot.title=element_text(hjust=0.5, vjust=0.9))
CorrMatrixPlot_PT1 <- CorrMatrixPlot

OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Row_Sums <- rowSums(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Row_Rank <-data.frame(rank(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Row_Sums))
colnames(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Row_Rank)<-"Rank"
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Rank <-cbind(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor, OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Row_Rank)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Rank_NCol <- ncol(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Rank)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order <- OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Rank[order(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Rank[,OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Rank_NCol]),]
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order$Rank <-NULL
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order <- OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order[,order(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order[nrow(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order),])]
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order<-as.matrix(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order <- melt(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order)
colnames(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order) <- c("X1", "X2", "Value")
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$X1 <- factor(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$X1, levels=unique(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$X1),ordered=TRUE)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$X2 <- factor(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$X2, levels=unique(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$X1),ordered=TRUE)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$Value<-cut(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$Value,breaks=c(-1,-0.60,-0.20,0.20,0.6,1),include.lowest=TRUE,label=c("-1.00 to -0.60 ","-0.60 to -0.20 ","-0.20 to 0.20 ","0.20 to 0.60 ","0.60 to 1.00 "))
CorrMatrix.NoPanel <- list(theme(panel.background=element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank())) 
CorrMatrixPlot_Order <- ggplot(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order,aes(x = X1, y = X2))
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
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot <-OAC_Converted_Transformed_Range_P_CORR
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number<-ncol(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot)
colnames(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot) <- c(1:OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor <- cor(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt <- melt(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor)
colnames(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt) <- c("X1", "X2", "Value")
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt$Value<-cut(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt$Value,breaks=c(-1,-0.60,-0.20,0.20,0.6,1),include.lowest=TRUE,label=c("-1.00 to -0.60 ","-0.60 to -0.20 ","-0.20 to 0.20 ","0.20 to 0.60 ","0.60 to 1.00 "))
CorrMatrix.NoPanel <- list(theme(panel.background=element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank())) 
CorrMatrixPlot <- ggplot(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt,aes(x = X1, y = X2))
CorrMatrixPlot <- CorrMatrixPlot + geom_tile(aes(fill=Value), colour="grey15")
CorrMatrixPlot <- CorrMatrixPlot + scale_fill_manual(values = c('-1.00 to -0.60 ' = "#B2182B", '-0.60 to -0.20 ' = "#F09B7A", '-0.20 to 0.20 ' = "#FFFFFF", '0.20 to 0.60 ' = "#87BEDA", '0.60 to 1.00 ' = "#2166AC"), name="Correlation ", drop=FALSE)
CorrMatrixPlot <- CorrMatrixPlot + scale_x_discrete(name="Variable Number", limits=c(1:OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number))
CorrMatrixPlot <- CorrMatrixPlot + scale_y_discrete(name="Variable Number", limits=c(1:OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number)) 
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

OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Row_Sums <- rowSums(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Row_Rank <-data.frame(rank(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Row_Sums))
colnames(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Row_Rank)<-"Rank"
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Rank <-cbind(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor, OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Row_Rank)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Rank_NCol <- ncol(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Rank)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order <- OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Rank[order(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Rank[,OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Rank_NCol]),]
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order$Rank <-NULL
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order <- OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order[,order(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order[nrow(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order),])]
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order<-as.matrix(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order <- melt(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order)
colnames(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order) <- c("X1", "X2", "Value")
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$X1 <- factor(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$X1, levels=unique(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$X1),ordered=TRUE)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$X2 <- factor(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$X2, levels=unique(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$X1),ordered=TRUE)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$Value<-cut(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$Value,breaks=c(-1,-0.60,-0.20,0.20,0.6,1),include.lowest=TRUE,label=c("-1.00 to -0.60 ","-0.60 to -0.20 ","-0.20 to 0.20 ","0.20 to 0.60 ","0.60 to 1.00 "))
CorrMatrix.NoPanel <- list(theme(panel.background=element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank())) 
CorrMatrixPlot_Order <- ggplot(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order,aes(x = X1, y = X2))
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
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot <-OAC_Converted_Transformed_Range_P_CORR
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number<-ncol(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot)
colnames(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot) <- c(1:OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor <- cor(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt <- melt(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor)
colnames(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt) <- c("X1", "X2", "Value")
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt$Value<-cut(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt$Value,breaks=c(-1,-0.60,-0.20,0.20,0.6,1),include.lowest=TRUE,label=c("-1.00 to -0.60 ","-0.60 to -0.20 ","-0.20 to 0.20 ","0.20 to 0.60 ","0.60 to 1.00 "))
CorrMatrix.NoPanel <- list(theme(panel.background=element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank())) 
CorrMatrixPlot <- ggplot(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt,aes(x = X1, y = X2))
CorrMatrixPlot <- CorrMatrixPlot + geom_tile(aes(fill=Value), colour="grey15")
CorrMatrixPlot <- CorrMatrixPlot + scale_fill_manual(values = c('-1.00 to -0.60 ' = "#B2182B", '-0.60 to -0.20 ' = "#F09B7A", '-0.20 to 0.20 ' = "#FFFFFF", '0.20 to 0.60 ' = "#87BEDA", '0.60 to 1.00 ' = "#2166AC"), name=" Correlation ", drop=FALSE)
CorrMatrixPlot <- CorrMatrixPlot + scale_x_discrete(name="Variable", limits=c(1:OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number),breaks = seq(1, OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number, 1), labels = c(1, rep("",OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number-2), OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number))
CorrMatrixPlot <- CorrMatrixPlot + scale_y_discrete(name="Variable", limits=c(1:OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number),breaks = seq(1, OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number, 1), labels = c(1, rep("",OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number-2), OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number))
CorrMatrixPlot <- CorrMatrixPlot + theme(axis.title.x = element_text(size=26,colour = "black", face = "bold"),axis.text.x=element_text(angle=0, size=20))
CorrMatrixPlot <- CorrMatrixPlot + theme(axis.title.y = element_text(size=26,colour = "black", face = "bold"),axis.text.y=element_text(size=20))+ CorrMatrix.NoPanel
CorrMatrixPlot <- CorrMatrixPlot + ggtitle("Correlation Matrix ") + theme(plot.title = element_text(size = 40, colour = "black", face = "bold"))
CorrMatrixPlot <- CorrMatrixPlot + coord_equal() + theme(legend.position = "bottom") 
CorrMatrixPlot <- CorrMatrixPlot + theme(legend.title = element_text(size=26)) + theme(legend.text = element_text(size=20)) 
CorrMatrixPlot <- CorrMatrixPlot + theme(axis.title.x=element_text(hjust=0.5, vjust=0.3))
CorrMatrixPlot <- CorrMatrixPlot + theme(axis.title.y=element_text(hjust=0.5, vjust=0.8))
CorrMatrixPlot <- CorrMatrixPlot + theme(plot.title=element_text(hjust=0.5, vjust=0.9))
CorrMatrixPlot_PT3 <- CorrMatrixPlot

OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Row_Sums <- rowSums(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Row_Rank <-data.frame(rank(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Row_Sums))
colnames(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Row_Rank)<-"Rank"
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Rank <-cbind(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor, OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Row_Rank)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Rank_NCol <- ncol(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Rank)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order <- OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Rank[order(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Rank[,OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Rank_NCol]),]
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order$Rank <-NULL
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order <- OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order[,order(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order[nrow(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order),])]
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order<-as.matrix(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order <- melt(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order)
colnames(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order) <- c("X1", "X2", "Value")
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$X1 <- factor(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$X1, levels=unique(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$X1),ordered=TRUE)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$X2 <- factor(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$X2, levels=unique(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$X1),ordered=TRUE)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$Value<-cut(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$Value,breaks=c(-1,-0.60,-0.20,0.20,0.6,1),include.lowest=TRUE,label=c("-1.00 to -0.60 ","-0.60 to -0.20 ","-0.20 to 0.20 ","0.20 to 0.60 ","0.60 to 1.00 "))
OAC_Converted_Transformed_Range_Cor_Order_Variables<-data.frame(unique(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$X1))
colnames(OAC_Converted_Transformed_Range_Cor_Order_Variables) <-"Ordered Variables"
OAC_Converted_Transformed_Range_Cor_Order_Variables_F<-OAC_Converted_Transformed_Range_Cor_Order_Variables[1,]
OAC_Converted_Transformed_Range_Cor_Order_Variables_L<-OAC_Converted_Transformed_Range_Cor_Order_Variables[ColNumber,]
OAC_Converted_Transformed_Range_Cor_Order_Variables_F <- as.matrix(OAC_Converted_Transformed_Range_Cor_Order_Variables_F)
OAC_Converted_Transformed_Range_Cor_Order_Variables_L <- as.matrix(OAC_Converted_Transformed_Range_Cor_Order_Variables_L)
OAC_Converted_Transformed_Range_Cor_Order_Variables_First <- as.numeric(OAC_Converted_Transformed_Range_Cor_Order_Variables_F)
OAC_Converted_Transformed_Range_Cor_Order_Variables_Last <- as.numeric(OAC_Converted_Transformed_Range_Cor_Order_Variables_L)
CorrMatrix.NoPanel <- list(theme(panel.background=element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank())) 
CorrMatrixPlot_Order <- ggplot(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order,aes(x = X1, y = X2))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + geom_tile(aes(fill=Value), colour="grey15")
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + scale_fill_manual(values = c('-1.00 to -0.60 ' = "#B2182B", '-0.60 to -0.20 ' = "#F09B7A", '-0.20 to 0.20 ' = "#FFFFFF", '0.20 to 0.60 ' = "#87BEDA", '0.60 to 1.00 ' = "#2166AC"), name=" Correlation ", drop=FALSE)
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + scale_x_discrete(name="Variable", labels = c(OAC_Converted_Transformed_Range_Cor_Order_Variables_First, rep("",OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number-2), OAC_Converted_Transformed_Range_Cor_Order_Variables_Last))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + scale_y_discrete(name="Variable",labels = c(OAC_Converted_Transformed_Range_Cor_Order_Variables_First, rep("",OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number-2), OAC_Converted_Transformed_Range_Cor_Order_Variables_Last))
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
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot <-OAC_Converted_Transformed_Range_P_CORR
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number<-ncol(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot)
colnames(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot) <- c(1:OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor <- cor(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt <- melt(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor)
colnames(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt) <- c("X1", "X2", "Value")
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt$Value<-cut(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt$Value,breaks=c(-1,-0.60,-0.20,0.20,0.6,1),include.lowest=TRUE,label=c("-1.00 to -0.60 ","-0.60 to -0.20 ","-0.20 to 0.20 ","0.20 to 0.60 ","0.60 to 1.00 "))
CorrMatrix.NoPanel <- list(theme(panel.background=element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank())) 
CorrMatrixPlot <- ggplot(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt,aes(x = X1, y = X2))
CorrMatrixPlot <- CorrMatrixPlot + geom_tile(aes(fill=Value))
CorrMatrixPlot <- CorrMatrixPlot + scale_fill_manual(values = c('-1.00 to -0.60 ' = "#B2182B", '-0.60 to -0.20 ' = "#F09B7A", '-0.20 to 0.20 ' = "#FFFFFF", '0.20 to 0.60 ' = "#87BEDA", '0.60 to 1.00 ' = "#2166AC"), name=" Correlation ", drop=FALSE)
CorrMatrixPlot <- CorrMatrixPlot + scale_x_discrete(name="Variable", limits=c(1:OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number),breaks = seq(1, OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number, 1), labels = c(1, rep("",OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number-2), OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number))
CorrMatrixPlot <- CorrMatrixPlot + scale_y_discrete(name="Variable", limits=c(1:OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number),breaks = seq(1, OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number, 1), labels = c(1, rep("",OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number-2), OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number))
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

OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Row_Sums <- rowSums(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Row_Rank <-data.frame(rank(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Row_Sums))
colnames(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Row_Rank)<-"Rank"
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Rank <-cbind(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor, OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Row_Rank)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Rank_NCol <- ncol(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Rank)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order <- OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Rank[order(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Rank[,OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Rank_NCol]),]
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order$Rank <-NULL
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order <- OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order[,order(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order[nrow(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order),])]
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order<-as.matrix(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order <- melt(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Cor_Order)
colnames(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order) <- c("X1", "X2", "Value")
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$X1 <- factor(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$X1, levels=unique(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$X1),ordered=TRUE)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$X2 <- factor(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$X2, levels=unique(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$X1),ordered=TRUE)
OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$Value<-cut(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$Value,breaks=c(-1,-0.60,-0.20,0.20,0.6,1),include.lowest=TRUE,label=c("-1.00 to -0.60 ","-0.60 to -0.20 ","-0.20 to 0.20 ","0.20 to 0.60 ","0.60 to 1.00 "))
OAC_Converted_Transformed_Range_Cor_Order_Variables<-data.frame(unique(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order$X1))
colnames(OAC_Converted_Transformed_Range_Cor_Order_Variables) <-"Ordered Variables"
OAC_Converted_Transformed_Range_Cor_Order_Variables_F<-OAC_Converted_Transformed_Range_Cor_Order_Variables[1,]
OAC_Converted_Transformed_Range_Cor_Order_Variables_L<-OAC_Converted_Transformed_Range_Cor_Order_Variables[ColNumber,]
OAC_Converted_Transformed_Range_Cor_Order_Variables_F <- as.matrix(OAC_Converted_Transformed_Range_Cor_Order_Variables_F)
OAC_Converted_Transformed_Range_Cor_Order_Variables_L <- as.matrix(OAC_Converted_Transformed_Range_Cor_Order_Variables_L)
OAC_Converted_Transformed_Range_Cor_Order_Variables_First <- as.numeric(OAC_Converted_Transformed_Range_Cor_Order_Variables_F)
OAC_Converted_Transformed_Range_Cor_Order_Variables_Last <- as.numeric(OAC_Converted_Transformed_Range_Cor_Order_Variables_L)
CorrMatrix.NoPanel <- list(theme(panel.background=element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank())) 
CorrMatrixPlot_Order <- ggplot(OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Melt_Order,aes(x = X1, y = X2))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + geom_tile(aes(fill=Value))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + scale_fill_manual(values = c('-1.00 to -0.60 ' = "#B2182B", '-0.60 to -0.20 ' = "#F09B7A", '-0.20 to 0.20 ' = "#FFFFFF", '0.20 to 0.60 ' = "#87BEDA", '0.60 to 1.00 ' = "#2166AC"), name=" Correlation ", drop=FALSE)
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + scale_x_discrete(name="Variable", labels = c(OAC_Converted_Transformed_Range_Cor_Order_Variables_First, rep("",OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number-2), OAC_Converted_Transformed_Range_Cor_Order_Variables_Last))
CorrMatrixPlot_Order <- CorrMatrixPlot_Order + scale_y_discrete(name="Variable",labels = c(OAC_Converted_Transformed_Range_Cor_Order_Variables_First, rep("",OAC_Converted_Transformed_Range_Correlation_Matrix_Plot_Var_Number-2), OAC_Converted_Transformed_Range_Cor_Order_Variables_Last))
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

CorrMatrixPlot_OrderLookup <- merge(OAC_Converted_Transformed_Range_Cor_Order_Variables, OAC_Converted_Transformed_Range_Input_Variables_Code, by.x = "Ordered Variables", by.y = "Variable Number", sort=F)
colnames(CorrMatrixPlot_OrderLookup) <-c("Variable Number", "Variable Code")
if(exists("OAC_Input_Lookup")=="TRUE")
{
OAC_Converted_Transformed_Range_Variable_Name_Lookup<-data.frame(OAC_Input_Lookup$VariableCode,OAC_Input_Lookup$VariableClustered)
CorrMatrixPlot_OrderLookup <- merge(CorrMatrixPlot_OrderLookup, OAC_Converted_Transformed_Range_Variable_Name_Lookup, by.x = "Variable Code", by.y = "OAC_Input_Lookup.VariableCode", sort=F)
CorrMatrixPlot_OrderLookup <- CorrMatrixPlot_OrderLookup[,c(2,1,3)]
colnames(CorrMatrixPlot_OrderLookup) <-c("Variable Number", "Variable Code","Variable Description")
}

OAC_Converted_Transformed_Range_Correlation_Matrix <-OAC_Converted_Transformed_Range_P_CORR
OAC_Converted_Transformed_Range_Correlation_Matrix_Var_Number<-ncol(OAC_Converted_Transformed_Range_Correlation_Matrix)
colnames(OAC_Converted_Transformed_Range_Correlation_Matrix) <- c(1:OAC_Converted_Transformed_Range_Correlation_Matrix_Var_Number)
OAC_Converted_Transformed_Range_Correlation_Matrix_Cor <- cor(OAC_Converted_Transformed_Range_Correlation_Matrix)
OAC_Converted_Transformed_Range_Correlation_Matrix_Melt <- melt(OAC_Converted_Transformed_Range_Correlation_Matrix_Cor)
OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_ID <- data.frame(1:nrow(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt))
colnames(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_ID) <-"ID"
OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Bind <-cbind(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_ID, OAC_Converted_Transformed_Range_Correlation_Matrix_Melt)
OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Abs <- data.frame(abs(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Bind))
OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr <- OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Abs[OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Abs$value <= 0.999999999999999,]
OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr <- OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr[OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr$value >=0.6,]
OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr <- OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr[with(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr, order(-value)), ]
OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr$value <- NULL
colnames(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr) <- c("ID", "X1", "X2")
colnames(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Bind)<- c("ID", "X1", "X2", "value")
OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr <- merge(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr, OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Bind, by= c("ID", "X1", "X2"), ,all = FALSE, sort=F)
OAC_Converted_Transformed_Range_Sig_Corr<- nrow(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr)
OAC_Converted_Transformed_Range_Sig_Corr_Range<-data.frame(1:OAC_Converted_Transformed_Range_Sig_Corr)
colnames(OAC_Converted_Transformed_Range_Sig_Corr_Range) <-("Order")
OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr <- cbind(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr, OAC_Converted_Transformed_Range_Sig_Corr_Range)
OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Plot<-OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr

Variable_Names<-data.frame(colnames(OAC_Converted_Transformed_Range_P_CORR))
Variable_Numbers<-data.frame(c(1:OAC_Converted_Transformed_Range_Correlation_Matrix_Var_Number))
Variable_LookupX1 <- cbind(Variable_Numbers, Variable_Names)
colnames(Variable_LookupX1) <- c("X1", "Variable 1")
Variable_LookupX2 <- cbind(Variable_Numbers, Variable_Names)
colnames(Variable_LookupX2) <- c("X2", "Variable 2")
Variable_Lookup <- merge(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr, Variable_LookupX1, by= "X1", ,all = FALSE, sort=F)
Variable_Lookup <- merge(Variable_Lookup, Variable_LookupX2, by= "X2", ,all = FALSE, sort=F)
Sig_Variable_Lookup<-Variable_Lookup[,c("Order","Variable 1", "Variable 2", "X1", "X2", "value")]
colnames(Sig_Variable_Lookup)<- c("ID", "Variable 1 Code", "Variable 2 Code", "Variable 1 Number", "Variable 2 Number", "Correlation Value")
Sig_Variable_Lookup_All<-Sig_Variable_Lookup[with(Sig_Variable_Lookup, order(ID)), ]
Sig_Variable_Lookup_Edit<-Sig_Variable_Lookup_All[Sig_Variable_Lookup_All$ID[Sig_Variable_Lookup_All$ID %% 2 ==0],]
Sig_Variable_Lookup_Edit_ID <-1:nrow(Sig_Variable_Lookup_Edit)
Sig_Variable_Lookup_Edit$ID <-NULL
Sig_Variable_Lookup_Edit<-cbind(Sig_Variable_Lookup_Edit_ID , Sig_Variable_Lookup_Edit)
colnames(Sig_Variable_Lookup_Edit)[1] <-"ID"

OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Plot$value<-cut(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr$value,breaks=c(-1,-0.9,-0.8,-0.7,0, 0.7,0.8,0.9,1),include.lowest=TRUE,label=c("-1.00 to -0.90 ","-0.90 to -0.80 ","-0.80 to -0.70 ","-0.70 to -0.60 ","0.60 to 0.70 ","0.70 to 0.80 ","0.80 to 0.90 ","0.90 to 1.00   "))

OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Var<-nrow(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr)
OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Var_X1 <-data.frame(unique(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Plot$X1))
OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Var_X2 <-data.frame(unique(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Plot$X2))
colnames(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Var_X1) <- ("UX1")
colnames(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Var_X2) <- ("UX2")
OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Var_X1 <- OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Var_X1[with(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Var_X1, order(UX1)), ]
OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Var_X2 <- OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Var_X2[with(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Var_X2, order(UX2)), ]
Sig_Corr_Var_X1_Min <- min(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Var_X1)
Sig_Corr_Var_X1_Max <- max(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Var_X1)
Sig_Corr_Var_X2_Min <- min(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Var_X2)
Sig_Corr_Var_X2_Max <- max(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Var_X2)
Sig_Corr_Var_Count<-count(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Var_X1)
Sig_Corr_Var_Count<-sum(Sig_Corr_Var_Count$freq)
Sig_Corr_Var_CNumX1<- data.frame(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Var_X1)
colnames(Sig_Corr_Var_CNumX1) <- ("X1")
Sig_Corr_Var_CNumX2<- data.frame(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Var_X2)
colnames(Sig_Corr_Var_CNumX2) <- ("X2")
Sig_Corr_Var_PNumZ1<- data.frame(1:Sig_Corr_Var_Count)
colnames(Sig_Corr_Var_PNumZ1) <- ("Z1")
Sig_Corr_Var_PNumZ2<- data.frame(1:Sig_Corr_Var_Count)
colnames(Sig_Corr_Var_PNumZ2) <- ("Z2")
Sig_Corr_VarX1 <-cbind(Sig_Corr_Var_CNumX1,Sig_Corr_Var_PNumZ1)
Sig_Corr_VarX2 <-cbind(Sig_Corr_Var_CNumX2,Sig_Corr_Var_PNumZ2)
OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Plot <- merge(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Plot, Sig_Corr_VarX1, by= "X1", ,all = FALSE, sort=F)
OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Plot <- merge(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Plot, Sig_Corr_VarX2, by= "X2", ,all = FALSE, sort=F)

CorrMatrixPlot_Sig <- ggplot(OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Plot,aes(x = Z1, y = Z2))
CorrMatrixPlot_Sig <- CorrMatrixPlot_Sig + geom_tile(aes(fill=value))
CorrMatrixPlot_Sig <- CorrMatrixPlot_Sig + scale_fill_manual(values = c('-1.00 to -0.90 ' = "#40004B", '-0.90 to -0.80 ' = "#762A83", '-0.80 to -0.70 ' = "#9970AB", '-0.70 to -0.60 ' = "#C2A5CF", '0.60 to 0.70 ' = "#A6DBA0", '0.70 to 0.80 ' = "#5AAE61", '0.80 to 0.90 ' = "#1B7837", '0.90 to 1.00   ' = "#00441B"), name="Correlation ", drop=FALSE)
CorrMatrixPlot_Sig <- CorrMatrixPlot_Sig + scale_x_discrete(name="Variable Number", limits=c(1:Sig_Corr_Var_Count),breaks = c(1:Sig_Corr_Var_Count), labels = OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Var_X1)
CorrMatrixPlot_Sig <- CorrMatrixPlot_Sig + scale_y_discrete(name="Variable Number", limits=c(1:Sig_Corr_Var_Count),breaks = c(1:Sig_Corr_Var_Count), labels = OAC_Converted_Transformed_Range_Correlation_Matrix_Melt_Sig_Corr_Var_X1)
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

dir.create("Correlation Data", showWarnings = FALSE)

write.table(OAC_Converted_Transformed_Range_Input_Variables, paste("Correlation Data/Correlation Matrix by Variables Order - ", OAC_Converted_Transformed_Range_Input_Variables[1,1], " to ", nrow(OAC_Converted_Transformed_Range_Input_Variables)," Variable Order.csv", sep = ""), sep = ",", row.names= FALSE, col.names = TRUE, qmethod = "double")

write.table(CorrMatrixPlot_OrderLookup, paste("Correlation Data/Correlation Matrix by Groupings Order - ", OAC_Converted_Transformed_Range_Cor_Order_Variables_First, " to ", OAC_Converted_Transformed_Range_Cor_Order_Variables_Last," Variable Order.csv", sep = ""), sep = ",", row.names= FALSE, col.names = TRUE, qmethod = "double")

CorrMatrixPlot_List<- list(CorrMatrixPlot_PT1, CorrMatrixPlot_PT2, CorrMatrixPlot_PT3, CorrMatrixPlot_PT4)
pdf(file = paste("Correlation Data/Correlation Matrix by Variables - ", ColNumber, " Variables - 4 Plot Types.pdf", sep=""), title = "Correlation Matrices by Variables", family='Courier', width=20, height=20)
CorrMatrixPlot_List
graphics.off()

CorrMatrixPlot_Order_List<- list(CorrMatrixPlot_Order_PT1, CorrMatrixPlot_Order_PT2, CorrMatrixPlot_Order_PT3, CorrMatrixPlot_Order_PT4)
pdf(file = paste("Correlation Data/Correlation Matrix by Groupings - ", ColNumber, " Variables - 4 Plot Types.pdf", sep=""), title = "Correlation Matrices by Groupings", family='Courier', width=20, height=20)
CorrMatrixPlot_Order_List
graphics.off()

CorrMatrixPlot_Sig_List<- list(CorrMatrixPlot_Sig_PT1, CorrMatrixPlot_Sig_PT2, CorrMatrixPlot_Sig_PT3)
pdf(file = paste("Correlation Data/Significant Correlation Matrix for ", ColNumber, " Variables - 3 Plot Types.pdf", sep=""), title = "Significant Correlation Matrices", family='Courier', width=20, height=20)
CorrMatrixPlot_Sig_List
graphics.off()

write.table(Sig_Variable_Lookup_All, paste("Correlation Data/Significant Correlation Matrix for ", ColNumber, " Variables Lookup Table - All.csv", sep = ""), sep = ",", row.names=FALSE, col.names = TRUE, qmethod = "double")

write.table(Sig_Variable_Lookup_Edit, paste("Correlation Data/Significant Correlation Matrix for ", ColNumber, " Variables Lookup Table - Edit.csv", sep = ""), sep = ",", row.names=FALSE, col.names = TRUE, qmethod = "double")

write.table(OAC_Converted_Transformed_Range_P_SEL_FINAL, paste("Correlation Data/Correlation and P Values - Pearson Correlation Coefficient - ", ColNumber, " Variables.csv", sep = ""), sep = ",", row.names=TRUE, col.names = NA, qmethod = "double")

write.table(OAC_Converted_Transformed_Range_S_SEL_FINAL, paste("Correlation Data/Correlation and P Values - Spearman Correlation Coefficient - ", ColNumber, " Variables.csv", sep = ""), sep = ",", row.names=TRUE, col.names = NA, qmethod = "double")

####################################################################################################
# Distribution Plots ###############################################################################
####################################################################################################

OAC_Converted_Transformed_Range_Plots <-data.frame(OAC_Converted_Transformed_Range_Input)
All_OA <-rownames(OAC_Converted_Transformed_Range_Plots)
YMin <- min(OAC_Converted_Transformed_Range_Plots)
YMax <- max(OAC_Converted_Transformed_Range_Plots)

NoR <-nrow(OAC_Converted_Transformed_Range_Plots)
NoC <-ncol(OAC_Converted_Transformed_Range_Plots)
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

OAC_Converted_Transformed_Range_PlotsSplit <- OAC_Converted_Transformed_Range_Plots
OA <-rownames(OAC_Converted_Transformed_Range_PlotsSplit)

for (ii in 1:PlotNoP)
{
Split<-get(paste("Split",ii, sep=""))
if (max(Split)> max(VarList))
{
SplitVarList <- min(Split):max(VarList)
assign(paste("OAC_Converted_Transformed_Range_Plots_",ii, sep=""), OAC_Converted_Transformed_Range_PlotsSplit[SplitVarList])
} else
{
assign(paste("OAC_Converted_Transformed_Range_Plots_",ii, sep=""), OAC_Converted_Transformed_Range_PlotsSplit[Split])
}
}

OAC_Converted_Transformed_Range_PlotsSplit <- cbind(OA,OAC_Converted_Transformed_Range_PlotsSplit)

for (v in 1:PlotNoP)
{
OAC_Converted_Transformed_Range_Plots_Selection<-get(paste("OAC_Converted_Transformed_Range_Plots_",v, sep=""))
assign(paste("OAC_Converted_Transformed_Range_Plots_",v, sep=""), cbind(OA,OAC_Converted_Transformed_Range_Plots_Selection))
}

OAC_Converted_Transformed_Range_PlotsSplit <- melt(OAC_Converted_Transformed_Range_PlotsSplit, id.vars="OA")
colnames(OAC_Converted_Transformed_Range_PlotsSplit) <-c("OA", "Variable", "Value")

for (g in 1:PlotNoP)
{
OACMeltSplit<-get(paste("OAC_Converted_Transformed_Range_Plots_",g, sep=""))
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
PlotSplit <- PlotSplit + geom_line(colour="black", size = PlotSplitMac, alpha = 0.8)
PlotSplit <- PlotSplit + scale_x_discrete(breaks=NULL) + theme_minimal(base_size = 12, base_family = "Courier") + theme(axis.ticks = element_blank(), axis.text.y = element_blank()) + theme(axis.title.y = element_blank(), axis.title.x = element_blank(),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank()) + scale_y_continuous(limits=c(YMin, YMax))
PlotSplit <- PlotSplit + facet_wrap( ~ Variable, ncol=FCol, nrow=FRow) + theme(strip.text.x = element_text(size=14))
assign(paste("PlotSplit_",h, sep=""), PlotSplit)
}

for (e in 1:PlotNoP)
{
OACMeltSplitPlot_Hist<-get(paste("OACMeltSplit_",e, sep=""))
PlotHist <- ggplot(OACMeltSplitPlot_Hist, aes(x=Value))
PlotHist<- PlotHist + geom_density() + theme_minimal(base_size = 12, base_family = "Courier") 
PlotHist<- PlotHist + theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(),axis.title.x = element_blank(),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
PlotHist<- PlotHist + facet_wrap(~ Variable, ncol=FCol, nrow=FRow, scales="free") + theme(strip.text.x = element_text(size=14))
assign(paste("PlotHist_",e, sep=""), PlotHist)
}

dir.create("Distribution Plots", showWarnings = FALSE)

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

OutputsEnd <- Sys.time()
StartSave <- Sys.time()

####################################################################################################
# Saving RData #####################################################################################
####################################################################################################

save.image("OAC_Converted_Transformed_Range_Clustered_with_Outputs.RData")

EndSave<- Sys.time()

OutputsEnd-OutputsStart
EndSave-StartSave

####################################################################################################
# End ##############################################################################################
####################################################################################################
