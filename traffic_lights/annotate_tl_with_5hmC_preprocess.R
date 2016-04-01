require(minfi)
raw_data <- read.table('GSE63179_series_matrix.txt', comment.char='!', header=TRUE, row.names=1)
grset <- makeGenomicRatioSetFromMatrix(as.matrix(raw_data))
gr<-getLocations(grset)
write.table(as.data.frame(gr), file = 'methregions_hg_19.txt')


