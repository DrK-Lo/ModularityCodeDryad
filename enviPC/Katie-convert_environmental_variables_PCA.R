
in1 <- read.table ("../enviPC/var_out_pine_all_COMBINED.table.contig_flt10.bayenv.envi2")
#?in1 <- read.table ("seqcap.bayenv.BF.envi2")
in2 <- t (in1)


nam1 <- read.table ("../enviPC/enviNamesAllAnalyses.txt")


in3 <- prcomp (in2, scale. = T)

(scores <- t (in3$x))


write.table (scores,"var_out_pine_all_COMBINED.table.contig_flt10.bayenv.envi2_PCA", col.names = F, row.names = F, quote = F)

write.table (scores,"seqcap.bayenv.BF.envi2_PCA", col.names = F, row.names = F, quote = F)