library(vegan)

data = read.table('views/projects/test/otus/denovo/table/table_filtered.csv', header=TRUE, sep='\t', row.names='OTUID')
meta = read.table('/bubo/home/h3/lucass/repos/illumitag/projects/test.csv', header=TRUE, sep='\t', row.names=1)

data_ordered = data[order(row.names(data)),]
data_sqrt = sqrt(data_ordered)
data_wa = wisconsin(data_sqrt)
data_dist = vegdist(data_wa, method="horn")

meta = meta[row.names(data),]
meta_ordered = meta[order(row.names(meta)),]
pool = factor(meta_ordered[,1])
chemistry = factor(meta_ordered[,3])

mod1 = betadisper(data_dist, chemistry)
permutest(mod1, control = permControl(nperm = 1000))
anova(mod1)
pdf(file="run_boxplot.pdf"); plot(mod1); boxplot(mod1); plot(TukeyHSD(mod1)); dev.off()

mod2 = betadisper(data_dist, pool)
permutest(mod2, control = permControl(nperm = 1000))
anova(mod2)
pdf(file="pool_boxplot.pdf"); plot(mod2); boxplot(mod2); plot(TukeyHSD(mod2)); dev.off()
