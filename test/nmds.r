library(vegan)
library(ggplot2)

data = read.table('views/projects/test/otus/denovo/table/table_filtered.csv', header=TRUE, sep='\t', row.names='OTUID')
meta = read.table('/bubo/home/h3/lucass/repos/illumitag/projects/test.csv', header=TRUE, sep='\t', row.names=1)

ord = metaMDS(data, distance='horn', trymax=200)


nmds = scores(ord)
df = merge(meta, nmds, by.x='row.names', by.y='row.names')
df$barcode = factor(df$barcode)
df$pool = factor(df$pool)
df$chemistry = factor(df$chemistry)

p = ggplot(df, aes(NMDS1, NMDS2)) + xlab("Dimension 1") + ylab("Dimension 2")
pdf(file="NMDS_pool.pdf"); p + geom_point(aes(colour=pool)); dev.off()
pdf(file="NMDS_barcode.pdf"); p + geom_point(aes(colour=barcode)); dev.off()
pdf(file="NMDS_chemistry.pdf"); p + geom_point(aes(colour=chemistry)); dev.off()
