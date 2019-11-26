library(dplyr)


#reading the txt file with classification
myCheckM <- read.table(file = "all_med_high_q_bins_checkM_stdout.txt", header = TRUE, 
                       col.names=c("Replicate", "Bin_ID", "Marker_Lineage", "UID",
                                   "genomes", "markers", "marker_sets", "_0", "_1", "_2", "_3",
                                   "_4", "_5", "Completeness",	"Contamination", "Strain_Heterogeneity"), sep= "")

#reading the csv file with sourmash results
myMash <- read.table(file = "SOUR_result.csv", header = TRUE, sep= ",")
#remove the number column, no meaning
myMash$X <- NULL

#reading the csv file with sourmash results
myANI <- read.table(file = "ANI_result.csv", header = TRUE, sep= ",")
#remove the number column, no meaning
myANI$X <- NULL

#creating the correct ID for matching
myCheckM$Replicate <- paste(myCheckM$Replicate, "_", sep="")
myCheckM$Bin_ID <- paste(myCheckM$Bin_ID, ".fa", sep="") 
myCheckM$ID <- do.call(paste, c(myCheckM[c("Replicate", "Bin_ID")], sep = "")) 

#ANI: creating classification for each cell, in a factor of 19
library(plyr)
myANI$CLASS_1 <- "N"
myANI$CLASS_2 <- "N"
myCheckM <- myCheckM[c("ID", "Marker_Lineage" )]
for (ad in myANI$AD_MAG_1){
  for (id in myCheckM$ID){
    if (ad == id){
      myANI[myANI$AD_MAG_1==ad,]$CLASS_1 <- myCheckM[myCheckM$ID==id,]$Marker_Lineage
    }
  }
}
for (ad in myANI$AD_MAG_2){
  for (id in myCheckM$ID){
    if (ad == id){
      myANI[myANI$AD_MAG_2==ad,]$CLASS_2 <- myCheckM[myCheckM$ID==id,]$Marker_Lineage
    }
  }
}

#ANI:Remove non-same lineage
myANI <- myANI[myANI$CLASS_1 == myANI$CLASS_2,]

#MASH: creating classification for each cell, in a factor of 19
myMash$CLASS_1 <- "N"
myMash$CLASS_2 <- "N"
for (ad in myMash$AD_MAG_1){
  for (id in myCheckM$ID){
    if (ad == id){
      myMash[myMash$AD_MAG_1==ad,]$CLASS_1 <- myCheckM[myCheckM$ID==id,]$Marker_Lineage
    }
  }
}
for (ad in myMash$AD_MAG_2){
  for (id in myCheckM$ID){
    if (ad == id){
      myMash[myMash$AD_MAG_2==ad,]$CLASS_2 <- myCheckM[myCheckM$ID==id,]$Marker_Lineage
    }
  }
}


#MASH:Remove non-same lineage
myMash <- myMash[myMash$CLASS_1 == myMash$CLASS_2,]


#Compute Unique Instances, account those in uniqueDiff in join
uniqueANI <- unique(myANI$AD_MAG_1)
uniqueMash <- unique(myMash$AD_MAG_1)
uniqueDiff1 <- setdiff (uniqueANI, uniqueMash)
uniqueDiff2 <- setdiff (uniqueMash, uniqueANI)

#Reduction
myANI <- myANI[c("AD_MAG_1", "AD_MAG_2")]
myMash <- myMash[c("AD_MAG_1", "AD_MAG_2")]

#Inner Join
join<-merge(x=myANI,y=myMash,by="AD_MAG_1")
colnames(join) <- c("AD_MAG_1", "AD_MAG_2", "AD_MAG_3")

#Changing factors to character
join$AD_MAG_1 <- as.character(join$AD_MAG_1)
join$AD_MAG_2 <- as.character(join$AD_MAG_2)
join$AD_MAG_3 <- as.character(join$AD_MAG_3)
join <- join[join$AD_MAG_2 == join$AD_MAG_3,]
join$AD_MAG_3 <- NULL

#Group into Directories
uniqueJoin <- unique(join$AD_MAG_1)
uniqueJoin <- as.character(uniqueJoin)
for(ad in uniqueJoin){
  join<-rbind(join, c(ad, ad))
}
join <- join[with(join,order(AD_MAG_1)),]


#Grouping
join$Group <- 0
count <- 1
first <- head(join,1)$AD_MAG_1
for(ad in join$AD_MAG_1){
  if (ad == first){
    join[join$AD_MAG_1==ad,]$Group <- count
  }
  else if (ad != first){
    first <- ad
    count <- count + 1
    join[join$AD_MAG_1==ad,]$Group <- count
  }
}
join$AD_MAG_1 <- NULL
colnames(join) <- c("MAG", "Group")

#Write to CSV
write.csv(join, "C:\\Users\\Eric\\Downloads\\consensus_bin\\consensus_bin_inner.csv", row.names=FALSE)

#Outer Join
join2<- rbind(myANI, myMash)
join2 <- join2[with(join2,order(AD_MAG_1)),]
join2 <- join2 %>% distinct()

#Grouping
join2$Group <- 0
count <- 1
first <- head(join2,1)$AD_MAG_1
for(ad in join2$AD_MAG_1){
  if (ad == first){
    join2[join2$AD_MAG_1==ad,]$Group <- count
  }
  else if (ad != first){
    first <- ad
    count <- count + 1
    join2[join2$AD_MAG_1==ad,]$Group <- count
  }
}
join2$AD_MAG_1 <- NULL
colnames(join2) <- c("MAG", "Group")

#Write to CSV
write.csv(join, "C:\\Users\\Eric\\Downloads\\consensus_bin\\consensus_bin_outer.csv", row.names=FALSE)

#Output in
#   \\consensus_bin_outer.csv: the result of outer join
#   \\consensus_bin_inner.csv: the result of inner join