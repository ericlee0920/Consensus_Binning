library(dplyr)
library(plyr)

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


#Outer Join
join2<- rbind(myANI, myMash)
join2 <- join2[with(join2,order(AD_MAG_1)),]
join2 <- join2 %>% distinct()

#Group into Directories
join2$AD_MAG_1 <- as.character(join2$AD_MAG_1)
join2$AD_MAG_2 <- as.character(join2$AD_MAG_2)
uniqueJoin <- unique(join2$AD_MAG_1)
uniqueJoin <- as.character(uniqueJoin)
for(ad in uniqueJoin){
  ad <- as.character(ad)
  join2<-rbind(join2, c(ad, ad))
}
join2 <- join2[with(join2,order(AD_MAG_1)),]


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



#ANI, MinHash: creating classification for each cell, in a factor of 19
join$CLASS <- "N"
join$CLASS_NAME <- "N"
join2$CLASS <- "N"
join2$CLASS_NAME <- "N"
myCheckM <- myCheckM[c("ID", "Marker_Lineage" )]
myCheckM$Marker_Lineage_2 <- as.character(myCheckM$Marker_Lineage)

for (ad in join$MAG){
  for (id in myCheckM$ID){
    if (ad == id){
      join[join$MAG==ad,]$CLASS <- myCheckM[myCheckM$ID==id,]$Marker_Lineage
      join[join$MAG==ad,]$CLASS_NAME <- myCheckM[myCheckM$ID==id,]$Marker_Lineage_2
      
    }
  }
}
for (ad in join2$MAG){
  for (id in myCheckM$ID){
    if (ad == id){
      join2[join2$MAG==ad,]$CLASS <- myCheckM[myCheckM$ID==id,]$Marker_Lineage
      join2[join2$MAG==ad,]$CLASS_NAME <- myCheckM[myCheckM$ID==id,]$Marker_Lineage_2
      
    }
  }
}


#create new line of median
freqfunc <- function(x, n){
  names(tail(sort(table(unlist(strsplit(as.character(x), ", ")))), n))
}
join$LINEAGE_MED <- 0
first <- head(join,1)$Group
for(ad in join$Group){
  d1<-join[join$Group == ad,]
  if (ad == first){
    join[join$Group == ad,]$LINEAGE_MED <- freqfunc(d1$CLASS, 1)
  }
  else if (ad != first){
    first <- ad
    join[join$Group == ad,]$LINEAGE_MED <- freqfunc(d1$CLASS, 1)
  }
}
join2$LINEAGE_MED <- 0
first <- head(join2,1)$Group
for(ad in join2$Group){
  d1<-join2[join2$Group == ad,]
  if (ad == first){
    join2[join2$Group == ad,]$LINEAGE_MED <- freqfunc(d1$CLASS, 1)
  }
  else if (ad != first){
    first <- ad
    join2[join2$Group == ad,]$LINEAGE_MED <- freqfunc(d1$CLASS, 1)
  }
}

#Write to CSV
write.csv(join, ".\\consensus_bin_inner.csv", row.names=FALSE)

#Write to CSV
write.csv(join2, ".\\consensus_bin_outer.csv", row.names=FALSE)



#CHECK: for all disagreement
count(join[join$CLASS!=join$LINEAGE_MED,])
freqfunc(join[join$CLASS!=join$LINEAGE_MED,]$CLASS, 3)
count(join2[join2$CLASS!=join2$LINEAGE_MED,])
freqfunc(join2[join2$CLASS!=join2$LINEAGE_MED,]$CLASS, 3)

#Joins: remove non-same lineage
join <- join[join$CLASS == join$LINEAGE_MED,]
join2 <- join2[join2$CLASS == join2$LINEAGE_MED,]


#Write to CSV
write.csv(join, ".\\consensus_bin_inner_conserve.csv", row.names=FALSE)

#Write to CSV
write.csv(join2, ".\\consensus_bin_outer_conserve.csv", row.names=FALSE)

#Output in
#   \\consensus_bin_outer.csv: the result of outer join
#   \\consensus_bin_inner.csv: the result of inner join
#   \\consensus_bin_outer_conserve.csv: the result of outer join with lineage trim
#   \\consensus_bin_inner_conserve.csv: the result of inner join with lingege trim
