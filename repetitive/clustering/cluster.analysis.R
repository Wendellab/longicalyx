library(ggplot2)
library(reshape2)
sessionInfo()

#################################
# R version 3.6.0 (2019-04-26)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 16299)
#
# Matrix products: default
#
# locale:
# [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
# [3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.1252    
#
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
#
# other attached packages:
# [1] reshape2_1.4.3 ggplot2_3.2.0 
#
# loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.1       crayon_1.3.4     withr_2.1.2      grid_3.6.0       plyr_1.8.4       gtable_0.3.0    
#  [7] magrittr_1.5     scales_1.0.0     pillar_1.4.1     stringi_1.4.3    rlang_0.4.0      lazyeval_0.2.2  
# [13] tools_3.6.0      stringr_1.4.0    munsell_0.5.0    compiler_3.6.0   pkgconfig_2.0.2  colorspace_1.4-1
# [19] tibble_2.1.3    
#################################

annot_clust <- read.table("AfricanRepeats.txt", header=T, sep="\t", row.names=1)


########### characterize composition ###########

# 8.5 multiplier represents # reads (x) * 85nt/read * 1 kb/1000nt * 100% = # reads * 8.5 = # Kb in entire genome for that class 
Kbamount <- data.frame(annot_clust[1], apply(annot_clust[2:ncol(annot_clust)], 2, function (x) x*8.5))
KBsum <- aggregate(. ~Lineage, data=Kbamount, FUN=sum)


for (species in unique(gsub("[_].*", "", names(KBsum[-1]))) ) {
    columns <- grep(species, names(KBsum))
    speciesSub <- subset(KBsum, select=as.numeric(columns))
    KBsum[,paste0(species,"mean")] <- rowMeans(speciesSub)
    KBsum[,paste0(species,"min")] <- apply(speciesSub,1,min)
    KBsum[,paste0(species,"max")] <- apply(speciesSub,1,max)
}

KBm <- melt(subset(KBsum, select=grep("mean|Lineage", names(KBsum))))
KBmin <- melt(subset(KBsum, select=grep("min", names(KBsum))))
KBmax <- melt(subset(KBsum, select=grep("max", names(KBsum))))

KBm$min <- KBmin$value
KBm$max <- KBmax$value

limits <- aes(ymax=KBm$max, ymin=KBm$min)

dodge <- position_dodge(width=0.9)

png("Figure_TE.amounts.png", 7500, 5000, pointsize=12, res=600)
ggplot(KBm, aes(x=Lineage, y=value, fill = variable)) + geom_bar(stat = "identity",position = dodge) + geom_errorbar(limits, position = dodge) + labs(title = "Aggregate amounts in each species", x="Broad element category", y="Aggregate amount (mean) in kilobases") + theme(axis.text = element_text(size = rel(1.5)), plot.margin=margin(2,2,2,2,"cm"), plot.title=element_text( face="bold", hjust=0.5), axis.title.x = element_text(face="bold", hjust=0.5), axis.title.y = element_text(face="bold", vjust=0.5))+theme_set(theme_grey(base_size=12))
dev.off()




