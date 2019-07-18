library(dplyr)

genomes <- read.table("F1.summary", sep = "\t", header=F, col.names=c("Genome","Chr","Family","Element","Length","Fragments","Copies","Solo_LTR","Total_Bp","Cover"))

genomes[is.na(genomes)] <- 0
TEsummary <- genomes %>%
    group_by(Genome,Family) %>%
    summarise(Fragments=sum(Fragments),Copies=sum(Copies),Total_Mb=sum(Total_Bp)/1000000,Coverage=mean(Cover))

TEdata <- as.data.frame(TEsummary %>% print(n = Inf))
TEdata <- TEdata[,-1]

write.table(TEdata, file = "TE.data", sep = "\t", col.names = T, quote = FALSE)
