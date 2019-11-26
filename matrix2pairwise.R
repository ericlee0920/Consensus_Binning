#OPTIONAL: install ggplot2
install.packages("ggplot2")

#run the pkg
library(ggplot2)

#reading the CSV file
myData <- read.table(file = "MAG_dist.csv", header = TRUE, row.names = 1, sep= ",")

#library for reshaping
library(reshape2)

#forming pairwise form
data <- melt(as.matrix(myData))

#renaming columns
names(data) <- c("AD_MAG_1", "AD_MAG_2", "JaccardIndex")

#renaming to avoid long names
data$AD_MAG_1<-substring(data$AD_MAG_1, 28)
data$AD_MAG_2<-substring(data$AD_MAG_2, 29)

#getting percentage needed for index, it's top 27612 samples that fit the 95 %
d1<-sum(data$JaccardIndex >= 0) * 0.006944605

#sort by JaccardIndex
data <- data[with(data,order(-JaccardIndex)),]

#take only top percentage
data <- data[1:27612,]

#removing Jaccard Index that is  0
data <- data[!(data$JaccardIndex == 0),]

#removing self comparisons
data<- data[data$AD_MAG_1!=data$AD_MAG_2,]

#writing to file
write.csv(data, "C:\\Users\\Eric\\Downloads\\sourmash\\data_original.csv")

#OPTIONAL: compute a histogram of `chol$AGE`
qplot(data$JaccardIndex, geom="histogram") 


#TODO:apply old Jaccard Index rule
#data <- data[data$JaccardIndex >= 0.65,]

#TODO:this is lowest similarity value in ANI
#41809	11A_II/BM_out.160.fa	AD152III_FD/BM_out.83.fa	95.0543
#3278122 AD152III_FD_BM_out.83.fa AD132III_FD_BM_out.106.fa which is only 0.156 in this


#reordering data
data <- data[order(data$AD_MAG_1),]

#copying data over
dataDeduplicate <-data

#deduplication, clouding
for (ad in dataDeduplicate$AD_MAG_1){
  d1<-dataDeduplicate[dataDeduplicate$AD_MAG_1==ad,]
  for (ad2 in d1$AD_MAG_2){
    dataDeduplicate<-dataDeduplicate[!(dataDeduplicate$AD_MAG_1==ad2),]
  }
}

#reordering data
dataDeduplicate <- dataDeduplicate[order(dataDeduplicate$AD_MAG_1),]

#remove
rm(d1, ad, ad2)

#writing to file
write.csv(dataDeduplicate, "C:\\Users\\Eric\\Downloads\\sourmash\\data_minhashed_dist.csv")

#removing myData
rm(myData, data)

