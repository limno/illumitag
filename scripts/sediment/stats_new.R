setwd("/Users/alper/Files/Illumina_tag/data")

library(vegan)
library(lmPerm)
library(ggplot2)
library(RSvgDevice)
library(gridExtra)
library(qvalue)


mdata<-read.csv(file="meta.data.csv", row.names=1, sep="\t")   ### sample data

d<-read.csv(file="otu_table.csv", sep="\t", row.names=1)


rd.m<-rrarefy(d,min(rowSums(d)))
dim(rd.m)
list2<-which(rowSums(d)>4000)
d2<-d[c(list2), ]
rd4000<-rrarefy(d2,min(rowSums(d2)))
dim(rd4000)

#nmds Soda Lakes
d1 = rd.m[ ,colSums(rd.m)!=0] 
dim(d1)
nmds1<-metaMDS(d1[c(1:9,34:42),],distance="bray", trymax=200)
plot(nmds1)

d.ill<-d1[c(1:10),] 
# just soda lakes d.ill<-d1[c(1:9),]
d.454<-d1[c(34:43),] 
# just soda lakes d.454<-d1[c(34:43),]

#procrustes test Soda Lakes
sqrt1<- sqrt(d.ill) 
w1<-wisconsin(sqrt1) 
a.ill<-vegdist(w1, method="bray")

sqrt2<- sqrt(d.454) 
w2<-wisconsin(sqrt2) 
a.454<-vegdist(w2, method="bray")

p<-procrustes(a.454,a.ill)
summary(p)
pt<-protest(a.454,a.ill)
mantel(a.454,a.ill)
plot(pt, type="t")
plot(a.ill,a.454)
a<-cbind(a.ill,a.454)

devSVG("regression_454_ill.svg")
ggplot(as.data.frame(a), aes(a.ill, a.454)) + geom_point(aes(size=1)) + theme_bw(10) + xlab("community dissimiarities based on pyrosequencing") + ylab("community dissimiarities based on illumina sequencing") + stat_smooth(method=lm, fullrange=TRUE) + geom_point()
dev.off()

am<-lm(a.ill~a.454)
summary(am)
pt

#wilcoxn test

dim(d.ill)
dim(d.454)

pvalues <- numeric(5111)
for (i in 1:5111) {
    pvalues[i] <- wilcox.test(x = d.ill[, i], y = d.454[, i],paired=TRUE, alternative="two.sided")$p.value
}

pvalues<-pvalues[!is.na(pvalues)]

qobj <- qvalue(pvalues)
q<-cbind(qobj$pvalue,qobj$qvalue)

l<-which(pvalues<0.05)
q[c(l),]

# nmds sediment
nmds2<-metaMDS(d1[c(10:33,43),],distance="bray", trymax=200)
plot(nmds2, type="t")


##### diversity estimates#################

eR<-as.data.frame(t(estimateR(rd4000)))  #chao1 and ACE
H <- diversity(rd4000)
S <- specnumber(rd4000)
J <- H/log(S)  #Pielou's
sh<-diversity(rd4000, index = "shannon", MARGIN = 1, base = exp(1)) #Shannon-Wiener
si<-diversity(rd4000, index = "simpson", MARGIN = 1, base = exp(1)) #Simpson

eR.ill<-eR[c(1:5,7,8),]
eR.454<-eR[c(32,33,35:38,40),]
chao1.lm<-lm(eR.ill$S.chao1~eR.454$S.chao1)
summary(chao1.lm)

devSVG("ill_454_chao1.svg")
plot(eR.454$S.chao1,eR.ill$S.chao1,cex=2)
abline(chao1.lm)
dev.off()


J.ill<-J[c(1:5,7,8)]
J.454<-J[c(32,33,35:38,40)]
J.lm<-lm(J.ill~J.454)
summary(J.lm)

devSVG("ill_454_J.svg")
plot(J.454,J.ill,cex=2)
abline(J.lm)
dev.off()


ACE.lm<-lm(eR.ill$S.ACE~eR.454$S.ACE)
summary(ACE.lm)

devSVG("ill_454_ACE.svg")
plot(eR.454$S.ACE,eR.ill$S.ACE,cex=2)
abline(ACE.lm)
dev.off()


sh.ill<-sh[c(1:5,7,8)]
sh.454<-sh[c(32,33,35:38,40)]
sh.lm<-lm(sh.ill~sh.454)
summary(sh.lm)

devSVG("ill_454_sh.svg")
plot(sh.454,sh.ill,cex=2)
abline(sh.lm)
dev.off()


si.ill<-si[c(1:5,7,8)]
si.454<-si[c(32,33,35:38,40)]
si.lm<-lm(si.ill~si.454)
summary(si.lm)

devSVG("ill_454_si.svg")
plot(si.454,si.ill,cex=2)
abline(si.lm)
dev.off()


### 250 sedminent samples ###

e<-read.csv(file="otu_table250.csv", sep="\t", row.names=1)
e4000<-e[rowSums(e) > 4000, ]
e4000<-rrarefy(e4000,min(rowSums(e4000)))
e4000 = e4000[ ,colSums(e4000)!=0] 

e5000<-e[rowSums(e) > 5000, ]
e5000<-rrarefy(e5000,min(rowSums(e5000)))
e5000 = e5000[ ,colSums(e5000)!=0] 

e1500<-e[rowSums(e) > 1500, ]
e1500<-rrarefy(e1500,min(rowSums(e1500)))
e1500 = e1500[ ,colSums(e1500)!=0] 

res4000<-merge(mdata,e4000,by.x="row.names", by.y="row.names")
adonis(formula = e4000 ~ pool * barcode, data = res4000, permutations = 1000, method = "bray")

res5000<-merge(mdata,e5000,by.x="row.names", by.y="row.names")
adonis(formula = e5000 ~ pool * barcode, data = res5000, permutations = 1000, method = "bray")

res1500<-merge(mdata,e1500,by.x="row.names", by.y="row.names")
adonis(formula = e1500 ~ as.factor(pool) * as.factor(barcode), data = res1500, permutations = 1000, method = "bray")

eR250<-as.data.frame(t(estimateR(e5000)))  #chao1 and ACE
res5000<-cbind(res5000,eR250)
H250 <- diversity(e5000)
S250 <- specnumber(e5000)
J250 <- H250/log(S250)  #Pielou's
res5000<-cbind(res5000,J250)
sh250<-diversity(e5000, index = "shannon", MARGIN = 1, base = exp(1)) #Shannon-Wiener
res5000<-cbind(res5000,sh250)
si250<-diversity(e5000, index = "simpson", MARGIN = 1, base = exp(1)) #Simpson
res5000<-cbind(res5000,si250)

anova.choa1<-aovp(S.chao1 ~ as.factor(pool) * as.factor(barcode),data=res5000)
anova.ACE<-aovp(S.ACE ~ as.factor(pool) * as.factor(barcode),data=res5000)
anova.J<-aovp(J250 ~ as.factor(pool) * as.factor(barcode),data=res5000)
anova.sh<-aovp(sh250 ~ as.factor(pool) * as.factor(barcode),data=res5000)
anova.si<-aovp(si250 ~ as.factor(pool) * as.factor(barcode),data=res5000)

summary(anova.choa1)
summary(anova.ACE)
summary(anova.J)
summary(anova.sh)
summary(anova.si)

an.chao<-aovp(S.chao1 ~ as.factor(pool),data=res5000)
TukeyHSD(an.chao)
an.ACE<-aovp(S.ACE ~ as.factor(pool),data=res5000)
TukeyHSD(an.ACE)
an.J<-aovp(J250 ~ as.factor(pool),data=res5000)
TukeyHSD(an.J)
an.sh<-aovp(sh250 ~ as.factor(pool),data=res5000)
TukeyHSD(an.sh)
an.si<-aovp(si250 ~ as.factor(pool),data=res5000)
TukeyHSD(an.si)

### pool 4 is different!
# comparing pools
e.1<-e1500[1:50,]
e.2<-e1500[51:100,]
e.3<-e1500[101:149,]
e.4<-e1500[150:198,]
e.5<-e1500[199:248,]
dist1<-metaMDSdist(e1500, distance="bray")
source<-c(rep(1,50),rep(2,50),rep(3,49),rep(4,49),rep(5,50))

mod1 <- betadisper(dist1, source)
mod1
permutest(mod1, control = permControl(nperm = 1000))
anova(mod1)
plot(mod1)

TukeyHSD(mod1)
plot(TukeyHSD(mod1))

boxplot(mod1)

adonis(formula = e.1 ~ c(1:50), permutations = 1000, method = "bray")
adonis(formula = e.2 ~ c(1:50), permutations = 1000, method = "bray")
adonis(formula = e.3 ~ c(1:49), permutations = 1000, method = "bray")
adonis(formula = e.4 ~ c(1:49), permutations = 1000, method = "bray")
adonis(formula = e.5 ~ c(1:50), permutations = 1000, method = "bray")

### clustering methods comparison

uparse<-read.csv(file="otu_table250.csv", sep="\t", row.names=1)
r.uparse<-rrarefy(uparse,min(rowSums(uparse)))
a.uparse<-metaMDSdist(r.uparse, distance="bray", trymax=200)

uclust<-read.csv(file="otu_table_uclust.csv", sep="\t", row.names=1)
r.uclust<-rrarefy(uclust,min(rowSums(uclust)))
a.uclust<-metaMDSdist(r.uclust, distance="bray", trymax=200)

cdhit<-read.csv(file="otu_table_cdhit.csv", sep="\t", row.names=1)
r.cdhit<-rrarefy(cdhit,min(rowSums(cdhit)))
a.cdhit<-metaMDSdist(r.cdhit, distance="bray", trymax=200)

t.test(rowSums(r.uparse>0),rowSums(r.uclust>0))
t.test(rowSums(r.uparse>0),rowSums(r.cdhit>0))
t.test(rowSums(r.uclust>0),rowSums(r.cdhit>0))

boxplot(rowSums(r.uparse>0),rowSums(r.cdhit>0),rowSums(r.uclust>0))

mean(rowSums(r.uparse>0))
mean(rowSums(r.uclust>0))
mean(rowSums(r.cdhit>0))

p1<-procrustes(a.uparse,a.uclust)
p2<-procrustes(a.uparse,a.cdhit)
p3<-procrustes(a.uclust,a.cdhit)

# absolute diversity
uparse1<-colSums(uparse)
uclust1<-colSums(uclust)
cdhit1<-colSums(cdhit)

uparse1.eR<-as.data.frame(estimateR(uparse1))
H <- diversity(uparse1)
S <- specnumber(uparse1)
uparse1.J <- H/log(S)  #Pielou's
uparse1.sh<-diversity(uparse1, index = "shannon", MARGIN = 1, base = exp(1)) #Shannon-Wiener
uparse1.si<-diversity(uparse1, index = "simpson", MARGIN = 1, base = exp(1)) #Simpson

source("rarefaction.txt")
uparse.rare<-rarefaction(uparse1, col=F)


uclust1.eR<-as.data.frame(estimateR(uclust1))
H <- diversity(uclust1)
S <- specnumber(uclust1)
uclust1.J <- H/log(S)  #Pielou's
uclust1.sh<-diversity(uclust1, index = "shannon", MARGIN = 1, base = exp(1)) #Shannon-Wiener
uclust1.si<-diversity(uclust1, index = "simpson", MARGIN = 1, base = exp(1)) #Simpson



#summary(p1)
pt1<-protest(a.uparse,a.uclust)
mantel(a.uparse,a.uclust)
plot(pt1, type="t")
plot(a.uparse,a.uclust)

#summary(p2)
pt2<-protest(a.uparse,a.cdhit)
mantel(a.uparse,a.cdhit)
plot(pt2, type="t")
plot(a.uparse,a.cdhit)

#summary(p3)
pt3<-protest(a.cdhit,a.uclust)
mantel(a.cdhit,a.uclust)
plot(pt3, type="t")
plot(a.cdhit,a.uclust)

# soda lakes comparison clustering methods

uparse1<-read.csv(file="otu_table_soda_uparse.csv", sep="\t", row.names=1)
r.uparse1<-rrarefy(uparse1,min(rowSums(uparse1)))
a.uparse1<-metaMDSdist(r.uparse1, distance="bray", trymax=200)

uclust1<-read.csv(file="otu_table_soda_uclust.csv", sep="\t", row.names=1)
r.uclust1<-rrarefy(uclust1,min(rowSums(uclust1)))
a.uclust1<-metaMDSdist(r.uclust1, distance="bray", trymax=200)

cdhit1<-read.csv(file="otu_table_soda_cdhit.csv", sep="\t", row.names=1)
r.cdhit1<-rrarefy(cdhit1,min(rowSums(cdhit1)))
a.cdhit1<-metaMDSdist(r.cdhit1, distance="bray", trymax=200)

t.test(rowSums(r.uparse1>0),rowSums(r.uclust1>0))
t.test(rowSums(r.uparse1>0),rowSums(r.cdhit1>0))
t.test(rowSums(r.uclust1>0),rowSums(r.cdhit1>0))

boxplot(rowSums(r.uparse1>0),rowSums(r.cdhit1>0),rowSums(r.uclust1>0))

mean(rowSums(r.uparse1>0))
mean(rowSums(r.uclust1>0))
mean(rowSums(r.cdhit1>0))


p4<-procrustes(a.uparse1,a.uclust1)
p5<-procrustes(a.uparse1,a.cdhit1)
p6<-procrustes(a.uclust1,a.cdhit1)

#summary(p1)
pt4<-protest(a.uparse1,a.uclust1)
mantel(a.uparse1,a.uclust1)
plot(pt4, type="t")
plot(a.uparse1,a.uclust1)

#summary(p2)
pt5<-protest(a.uparse1,a.cdhit1)
mantel(a.uparse1,a.cdhit1)
plot(pt5, type="t")
plot(a.uparse1,a.cdhit1)

#summary(p3)
pt6<-protest(a.cdhit1,a.uclust1)
mantel(a.cdhit1,a.uclust1)
plot(pt6, type="t")
plot(a.cdhit1,a.uclust1)

uparse5000<-uparse1[rowSums(uparse1) > 5000, ]
uparse5000<-rrarefy(uparse5000,min(rowSums(uparse5000)))
uparse5000 = uparse5000[ ,colSums(uparse5000)!=0] 

uclust5000<-uclust1[rowSums(uclust1) > 5000, ]
uclust5000<-rrarefy(uclust5000,min(rowSums(uclust5000)))
uclust5000 = uclust5000[ ,colSums(uclust5000)!=0] 

cdhit5000<-cdhit1[rowSums(cdhit1) > 5000, ]
uclust5000<-rrarefy(uclust5000,min(rowSums(uclust5000)))
uclust5000 = uclust5000[ ,colSums(uclust5000)!=0] 

eRuparse<-as.data.frame(t(estimateR(uparse5000)))  #chao1 and ACE
Huparse <- diversity(uparse5000)
Suparse <- specnumber(uparse5000)
Juparse <- Huparse/log(Suparse)  #Pielou's

eRuclust<-as.data.frame(t(estimateR(uclust5000)))  #chao1 and ACE
Huclust <- diversity(uclust5000)
Suclust <- specnumber(uclust5000)
Juclust <- Huclust/log(Suclust)  #Pielou's

eRcdhit<-as.data.frame(t(estimateR(cdhit5000)))  #chao1 and ACE
Hcdhit <- diversity(cdhit5000)
Scdhit <- specnumber(cdhit5000)
Jcdhit <- Hcdhit/log(Scdhit)  #Pielou'

eRsoda<-merge(eRuparse,eRuclust,by.x="row.names", by.y="row.names")
eR.pc<-lm(eRsoda$S.chao1.x~eRsoda$S.chao1.y)
summary(eR.pc)

eRsoda<-merge(eRuparse,eRcdhit,by.x="row.names", by.y="row.names")
eR.ph<-lm(eRsoda$S.chao1.x~eRsoda$S.chao1.y)
summary(eR.ph)

eRsoda<-merge(eRcdhit,eRuclust,by.x="row.names", by.y="row.names")
eR.hc<-lm(eRsoda$S.chao1.x~eRsoda$S.chao1.y)
summary(eR.hc)

Jsoda<-merge(as.matrix(Juparse), as.matrix(Juclust),by.x="row.names", by.y="row.names")
J.pc<-lm(Jsoda$V1.x~Jsoda$V1.y)
summary(J.pc)

Jsoda<-merge(as.matrix(Juparse), as.matrix(Jcdhit),by.x="row.names", by.y="row.names")
J.ph<-lm(Jsoda$V1.x~Jsoda$V1.y)
summary(J.ph)

Jsoda<-merge(as.matrix(Jcdhit), as.matrix(Juclust),by.x="row.names", by.y="row.names")
J.hc<-lm(Jsoda$V1.x~Jsoda$V1.y)
summary(J.hc)

# betadispersion pools
d<-rep(1,248)
mod1 <- betadisper(a.uparse,d, type="centroid")
mod2 <- betadisper(a.cdhit,d, type="centroid")
mod3 <- betadisper(a.uclust,d, type="centroid")

boxplot(mod1$distance, mod2$distance, mod3$distance)

# mean distances of uparse matrix
dd1<-metaMDSdist(d1,distance="bray", trymax=200)
dd1<-as.matrix(dd1)
mean(dd1[10:33,10:33]) # MEAN of replicated sediment
mean(dd1[10:33,43]) # MEAN of 454 to Illumina replicates
dd2<-melt(dd1)[melt(upper.tri(dd1))$value,]
dd2<-dd2[order(dd2$Var1),]
row.names(dd2)<-c(1:903)

mean(d[c(9:32,43:49,84:89,124:128,163:166,201:203,238,239,274),3])

mean(dd2[c(9:32,50:73,90:113,129:152,167:190,204:227,240:263,275:298,309:332),3])





#p3<-ggplot(a, aes(River_km, J)) + geom_point(aes(size=1,colour=as.factor(Filter_fraction))) + theme_bw(10) + xlab("River kilometer") + ylab("J") + scale_fill_discrete(name="Filter\nsize") + stat_smooth(method=lm, fullrange=TRUE, aes(fill = factor(Filter_fraction))) + geom_point()

#p4<-ggplot(b, aes(River_km, eR$S.chao1)) + geom_point(aes(size=1,colour=as.factor(Filter_fraction))) + theme_bw(10) + xlab("River kilometer") + ylab("chao1") + scale_fill_discrete(name="Filter\nsize") + stat_smooth(method=lm, fullrange=TRUE, aes(fill = factor(Filter_fraction))) + geom_point()


