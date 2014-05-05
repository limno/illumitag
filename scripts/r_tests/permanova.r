library(vegan)

data = read.table('views/projects/test/otus/denovo/table/table_filtered.csv', header=TRUE, sep='\t', row.names='OTUID')
meta = read.table('/bubo/home/h3/lucass/repos/illumitag/projects/test.csv', header=TRUE, sep='\t', row.names=1)

data_ordered = data[order(row.names(data)),]

meta = meta[row.names(data),]
meta_ordered = meta[order(row.names(meta)),]
meta_ordered$pool = factor(meta_ordered$pool)
meta_ordered$barcode = factor(meta_ordered$barcode)
meta_ordered$chemistry = factor(meta_ordered$chemistry)
adonis(formula = data_ordered ~ pool * barcode, data=meta_ordered, permutations=1000, method="horn")

subset = row.names(meta_ordered[meta_ordered$pool == 1 | meta_ordered$pool == 3,])
subdata = data_ordered[subset,]
submeta = meta_ordered[subset,]
adonis(formula = subdata ~ pool * barcode, data=submeta, permutations=1000, method="horn")