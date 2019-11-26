#reading the tab delimited file
data <- read.table("results.txt", header = TRUE, sep= "")

#removing ref and hit columns
data$Ref <- NULL
data$Hit <- NULL

#counting instances above 95%, it's 0.006944605 of total samples
sum(data$Similarity >= 95) / (1994*1994)

#applying 95% ANI similarity rule
data <- data[data$Similarity >= 95,]

#reordering data
data <- data[order(data$AD_MAG_1, data$AD_MAG_2),]

#OPTIONAL: forming matrix
dfr <- reshape(data, direction = "wide", timevar = "AD_MAG_1", idvar = "AD_MAG_2")

#OPTIONAL: setting na to 0
dfr[is.na(dfr)] <- 0

#OPTIONAL: writing matrix to file [THIS MATRIX IS NOT DEDUPLICATED]
write.csv(dfr, "C:\\Users\\Eric\\Downloads\\fastANI\\data_95.csv")

#removing self comparisons
dataDeduplicate<- data[data$AD_MAG_1!=data$AD_MAG_2,]

#deduplication, clouding
for (ad in dataDeduplicate$AD_MAG_1){
  d1<-dataDeduplicate[dataDeduplicate$AD_MAG_1==ad,]
  for (ad2 in d1$AD_MAG_2){
    dataDeduplicate<-dataDeduplicate[!(dataDeduplicate$AD_MAG_1==ad2),]
  }
}

#renaming to avoid long names
dataDeduplicate$AD_MAG_1<-substring(dataDeduplicate$AD_MAG_1, 28)
dataDeduplicate$AD_MAG_2<-substring(dataDeduplicate$AD_MAG_2, 28)

#writing to file
write.csv(dataDeduplicate, "C:\\Users\\Eric\\Downloads\\fastANI\\data_deduplicated.csv")

#removing myData
rm(myData)
rm(ad, ad2, d1, data)
