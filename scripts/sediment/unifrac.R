setwd("/Users/alper/Files/Illumina_tag/data")

library(vegan)
library(lmPerm)
library(ggplot2)
library(RSvgDevice)
library(gridExtra)
library(reshape2)


mdata<-read.csv(file="meta.data.csv", row.names=1, sep="\t")   ### sample data

eval<-read.csv(file="distance_evaluation.csv", sep="\t", row.names=1)
pyro<-read.csv(file="distance_pyro.csv", sep="\t", row.names=1)
soda<-read.csv(file="distance_soda.csv", sep="\t", row.names=1)



nmds.eval<-metaMDS(eval, trymax=200)
plot(nmds.eval)

nmds.pyro<-metaMDS(pyro, trymax=200)
plot(nmds.pyro, type="t")

nmds.soda<-metaMDS(soda, trymax=200)
plot(nmds.pyro)


### 250 sedminent samples ###


res<-merge(mdata,eval,by.x="row.names", by.y="row.names")
adonis(formula = eval ~ pool * barcode, data = res, permutations = 1000)

### comparison 454 vs. Ilumina

rownames(pyro)
d.ill<-pyro[c(1:10),] 
d.454<-pyro[c(34:43),] 

p<-procrustes(as.dist(d.454),as.dist(d.ill))
summary(p)
pt<-protest(as.dist(d.454),as.dist(d.ill))
mantel(as.dist(d.454),as.dist(d.ill))
plot(pt, type="t")
plot(as.dist(d.454),as.dist(d.ill))
a<-as.data.frame(cbind(as.dist(d.454),as.dist(d.ill)))

devSVG("regression_454_ill_unifrac.svg")
ggplot(a, aes(V1,V2)) + geom_point(aes(size=1)) + theme_bw(10) + xlab("unifrac distance based on pyrosequencing") + ylab("unifrac distance based on illumina sequencing") + stat_smooth(method=lm, fullrange=TRUE) + geom_point()
dev.off()

# mean distances

mean(as.matrix(pyro[10:33,10:33])) # MEAN of replicated sediment
mean(as.matrix(pyro[10:33,43])) # MEAN of 454 to Illumina replicates
d<-melt(as.matrix(pyro))[melt(upper.tri(as.matrix(pyro)))$value,]
d<-d[order(d$Var1),]
row.names(d)<-c(1:903)

mean(d[c(9:32,43:49,84:89,124:128,163:166,201:203,238,239,274),3])

mean(d[c(1:8,50:73,90:113,129:152,167:190,204:227,240:263,275:298,309:332),3])
