
rm(list = ls())

library("ggplot2")
library("reshape2")
library("cowplot")
library("wesanderson")
library("png")
library("grid")

se <- function(x){sd(x)/sqrt(length(x))}

pal <- wes_palette("Zissou1", 3, type = "continuous")
pal2 <- wes_palette("Moonrise2", 4, type = "continuous")

source("R/data_prep.R")

dat <- read.csv("data/traits.csv")
spp <- read.csv("data/species.csv")
env <- read.csv("data/envs.csv")
lit <- read.csv("data/abundance.csv")

dat$Location <- factor(dat$Location, levels=c("Outer","Lagoon","Inner"))

dat$spp <- spp$spp[match(dat$Species, spp$species)]
dat$Occ <- spp$Occ[match(dat$Species, spp$species)]
dat$SiteID <- env$SiteID[match(dat$Site, env$Site)]
dat$Occ2 <- ifelse(dat$Occ=="Both", "Generalists", "Specialists")

# shades
shades <- dat$col
names(shades) <- dat$ID


####################################
####################################
####################################
# ENVIRONMENTS

envlong <- melt(subset(env, select=-c(Date, X, Depth_m, long, long2, lat, SpCond_mS.cm)), id.var=c("Site", "Location", "SiteID"))
head(envlong)
 
 ggplot(envlong, aes(Location, value))+
 geom_boxplot(aes(fill=Location))+
 facet_wrap(~variable, scales="free_y")

####################################
####################################
####################################
# ABUNDANCE

nrow(lit)

a_tran <- lit
a_tran <- aggregate(.~SiteID+Location+Site+Circle+ Transect+Species, a_tran, sum) 

head(a_tran)
a_site <- aggregate(.~SiteID+Species+Site+Location, subset(lit, select=-c(Circle, Transect, X)), mean)

UNSAMP<-as.character(unique(a_site$Species[!a_site$Species %in% c(spp$species, "Porites Massive")]))

a_site$SpeciesS <- ifelse(a_site$Species %in% UNSAMP, "Other", as.character(a_site$Species))

ggplot(a_site, aes(x=SiteID, y=cover, group=Species, 
fill=SpeciesS))+ 
geom_area(position = 'stack', col="black", size=0.1)

####################################
####################################
####################################
# OCCURENCE

a_site$Occ <- spp$Occ[match(a_site$Species, spp$species)]

ggplot()+
geom_path(data=a_site[a_site$cover > 0,], aes(SiteID, SpeciesS, group=SpeciesS))+guides(col="none")+
geom_point(data=dat, aes(SiteID, Species, col=Species))+
facet_grid(Occ~., scales="free_y", space="free_y")

####################################
####################################
####################################
#  % accounted for in sampling
# linking groups

a_site$samp <- ifelse(a_site$SpeciesS=="Unsampled", "n", "y")
#a_site$samp <- ifelse(a_site$SpeciesS=="Porites Massive", "n", a_site$samp)
total <- aggregate(cover~SiteID, a_site, sum)
sampd <- aggregate(cover~samp+SiteID, a_site, sum)
sampd$tot <- total$cover[match(sampd$SiteID, total$SiteID)]
sampd$prop <- sampd$cover/sampd$tot *100
sampd
aggregate(prop~samp, sampd, mean)
 
####################################
####################################
####################################
# TRAIT PCA

traits <- c("chl_a","zoox", "protein", "SD_gcm3", "branch_density","biomass", "SAV", "polyp_size",  "branch_width")

info <- read.csv("data/raw_data/trait_info.csv")

tdat <- dat 
tdat$Occ<-as.character(tdat$Occ)
tdat <- tdat[tdat$Occ=="Both",]

pca <- prcomp(na.omit(log(tdat[,traits])), center=T, scale=T)
biplot(pca)

####################################
####################################
####################################
# Fig.1

source("R/map_palau.R")
source("figs/fig1.R")
fig.1
#ggsave("figs/fig.1.jpg", fig.1, width=140, height=108, units="mm", dpi = 1000)

####################################
####################################
####################################
# effect sizes .... 

library("lme4")
library("MuMIn")
library("lmerTest")

cohensD <- function(B){
	# lme4 cohens d=(2*t)/sqrt(df), where t = b/se(b)
	tval <- B/coef(summary(mod))["x","Std. Error"]
	2*tval/sqrt(coef(summary(mod))["x","df"])	
	}

params <- NULL
for (t in traits){
	tdat$y <- scale(log(tdat[,t]))
	tdat$x <- ifelse(tdat$Location=="Inner", 1,ifelse(tdat$Location=="Outer",0, NA))
	df <- tdat[!is.na(tdat$x),]
	df <- df[!is.na(df$y),]
	mod <- lmer(y ~ x + (1|Species), data=df)
	coefs <- fixef(mod)
	dscore <- cohensD(coefs[2])
	# EMAtools::lme.dscore(mod, df, "lme4") #check
	CI.x <- confint(mod, method="Wald")
	d.upper <-cohensD(CI.x ["x","97.5 %"])
	d.lower <-cohensD(CI.x ["x","2.5 %"])
	rsq.m <- r.squaredGLMM(mod)[1]
	rsq.c <- r.squaredGLMM(mod)[2]
	co <- coef(summary(mod))
	se.x <- co["x","Std. Error"]
	pval.x <- co["x","Pr(>|t|)"]
	T.x <- co["x","t value"]
	df.x <- co["x","df"]
	int <- co["(Intercept)","Estimate"]
	se.i <- co["(Intercept)","Std. Error"]
	pval.i <-co["(Intercept)","Pr(>|t|)"]
	T.i <- co["(Intercept)","t value"]
	df.i <- co["(Intercept)","df"]
	params <- rbind(params, data.frame(t=t, x=coefs[2], se.x, pval.x, int, pval.i, se.i, T.x, T.i,df.i, df.x, rsq.m, rsq.c,  sd.resid=sd(residuals(mod, type="response")),  row.names=NULL, dscore, d.upper, d.lower))
	}
		
params$e.size <- params$x / params$sd.resid
params$e.se <- params$se.x / params$sd.resid *1.96 
params$sig <- ifelse(params$pval.x < 0.05, '*', "")
params$label<- info$label[match(params$t, info$trait)]
params$type <- info$type[match(params$t, info$trait)]
params$name <- info$name[match(params$t, info$trait)]

ggplot(data=params)+
geom_vline(xintercept=0, col="grey")+
geom_point(aes(x=e.size, y=reorder(label, e.size)))+
geom_segment(aes(x=e.size+e.se, xend=e.size-e.se, y=reorder(label, e.size), yend=reorder(label, e.size)))+
geom_point(aes(x=dscore, y=reorder(label, e.size)), col="red", , size=1)+
geom_segment(aes(x=d.upper, xend=d.lower, y=reorder(label, e.size), yend=reorder(label, e.size)), col="red", size=0.5)
 
####################################
####################################
####################################
 
# species specific responses
# change in averages inside to out

splist <- unique(tdat$Species[tdat$Occ=="Both"])

# split lagoon data for sample sizes
tdat2 <- dat
tdat2$Location <- ifelse(tdat2$Species=="Acropora cf. subglabra" & tdat2$SiteID=="E", "Outer", as.character(tdat2$Location))
tdat2$Location <- ifelse(tdat2$Species=="Acropora cf. muricata" & tdat2$SiteID=="D", "Inner", as.character(tdat2$Location))

params2 <- NULL

for(sp in splist){
	for(t in traits){
df <- tdat2[tdat2$Species==sp,]
df$y <- df[,t]
df$x <- ifelse(df$Location=="Inner", 1, ifelse(df$Location=="Outer",0, NA))
mod <- lm(y ~ x, df)
summary(mod)
coefs <- coef(mod)
se.x <- coef(summary(mod))["x","Std. Error"]
pval.x <- coef(summary(mod))["x","Pr(>|t|)"]
sig <- ifelse(pval.x<0.05, "*","")
params2 <- rbind(params2, data.frame(x=coefs[2], se.x, pval.x, sd.resid=sd(residuals(mod, type="response")), t=t, sp=sp, sig, row.names=NULL))
   }
  }

params2$esize <- params2$x/params2$sd.resid
params2$type <- info$type[match(params2$t, info$trait)]
params2$label <- info$label[match(params2$t, info$trait)]
params2$name <- info$name[match(params2$t, info$trait)]
params2$spp<- spp$spp[match(params2$sp, spp$species)]
params2$sig2 <- ifelse(params2$sig=="", "ns", as.character(params2$type))
params2$e.se <- params2$se.x / params2$sd.resid *1.96

params2$label <- factor(params2$label, levels=c("ZD","CH","PB", "TB", "SD", "SV", "BD", "CW", "BW"))

ggplot(params2, aes(label, esize))+
geom_bar( stat="identity", aes(fill=type))+
facet_wrap(~spp,scale="free_x")+
coord_flip()+geom_hline(yintercept=0, size=0.1)+
lims(y=c(-4.5,4))+
labs(y="Effect size", x="Trait")

####################################
####################################
####################################
# Figure 2...
source("figs/fig2.R")
fig.2
#ggsave("figs/fig.2.jpeg", fig.2, width=215, height=100, units="mm", dpi = 1000)

####################################
####################################
####################################
# Locally dominant species 

dat$spp <- spp$spp[match(dat$Species, spp$species)]
pca3 <- prcomp(na.omit(log(dat[,traits])), center=T, scale=T)
biplot(pca3)

pcdat <- cbind(na.omit(dat[,c(traits, "ID","Species","Species_TB", "Genus","SiteID","Location", "Occ", "spp")]), data.frame(pc1=pca3$x[,1], pc2=pca3$x[,2]))

####################################
####################################
####################################

# CWMs
library("FD")

# SIMPLEST FORM OF CWMS - global mean
dat2 <- dat
dat2$pc1<-pcdat$pc1[match(dat2$ID, pcdat$ID)]
dat2$pc2<-pcdat$pc2[match(dat2$ID, pcdat$ID)]
dat2$Species[dat2$Species %in% c("Lobophyllia cf. hemprichii", "Lobophyllia cf. corymbosa")]<- "Lobophyllia spp."
dat2$Species[dat2$Species %in% c("Porites cf. attenuata")]<- "Porites cf. cylindrica"
nrow(dat2)

fixed <- aggregate(.~Species, dat2[,c("Species", traits, "pc1", "pc2")], mean, na.rm=T, na.action=NULL)
fixed <- data.frame(fixed[,c(traits,"pc1","pc2")], row.names=fixed$Species)
head(fixed)

# Abundances...
abun <- a_tran
abun <- abun[abun$Species %in% unique(dat2$Species),]
abun$ID <- paste(abun$SiteID, abun$Circle, abun$Transect, sep=".")
abun2 <- acast(abun, ID~Species, value.var="cover")
abun2 <- abun2[,rownames(fixed)]
head(abun2)
rownames(fixed)
colnames(abun2)==rownames(fixed)

f.cwm <- functcomp(a=as.matrix(abun2), x=fixed)
f.cwm$FDis <- fdisp(a=as.matrix(abun2), d=dist(fixed))$FDis

####################################
####################################
####################################

# CWMS with site variations... 
specific <- aggregate(.~Species+SiteID, dat2[,c("Species","SiteID", traits, "pc1","pc2")], function(x){mean(x, na.rm=T)}, na.action=na.pass)
head(specific)

specific_loc <- aggregate(.~Species+Location, dat2[,c("Species","Location", traits, "pc1","pc2")], function(x){mean(x, na.rm=T)}, na.action=na.pass)
head(specific_loc)

s.cwm <- NULL
for (site in unique(dat2$SiteID)){
loc <- dat$Location[dat2$SiteID==site][1]
# abundances
site.a <- subset(abun, SiteID==site)
site.a <- acast(site.a, ID~Species, value.var="cover")
site.a <- site.a[,!colSums(site.a)==0] # abundances
# traits
site.x <- subset(specific, SiteID==site)
site.x <- data.frame(site.x[,c(traits,"pc1","pc2")], row.names=site.x$Species)
# infill traits by location if no sample
got <- rownames(site.x)[rownames(site.x) %in% colnames(site.a)]
notgot <- colnames(site.a)[!colnames(site.a) %in% got]
site.x1 <- site.x[got,]
site.x2 <- subset(specific_loc, Location==loc)
site.x2 <- data.frame(site.x2[,c(traits,"pc1","pc2")], row.names=site.x2$Species)
site.x3 <- site.x2[notgot,]
site.x4 <- rbind(site.x1, site.x3)[colnames(site.a),]
fd <- fdisp(d=dist(site.x4), a=as.matrix(site.a))$FDis
s.cwm <- rbind(s.cwm,  cbind(functcomp(a=as.matrix(site.a), x=site.x4), FDis=fd))
}

head(s.cwm)
head(f.cwm)

####################################
# CWM change

cdat <- rbind(cbind(s.cwm, cwm="specific"), cbind(f.cwm, cwm="fixed"))
cdat$ID <- substr(rownames(cdat), 1,5)
cdat$Location <- abun$Location[match(cdat$ID, abun$ID)] 
cdat$SiteID <- abun$SiteID[match(cdat$ID, abun$ID)] 
head(cdat)

cdat2 <- melt(cdat, variable.name="trait", value.name="y")
cdat2 <- cdat2[cdat2$cwm=="specific",]
maxvals <- aggregate(y~trait, cdat2, max)
cdat2$max <- maxvals$y[match(cdat2$trait, maxvals$trait)]
cdat2$y_norm <- cdat2$y/cdat2$max
cdat2$type <- info$type[match(cdat2$trait, info$trait)]
cdat2$label <- info$label[match(cdat2$trait, info$trait)]
head(cdat2)

cdat_av <- aggregate(y_norm~SiteID+trait+type+label, cdat2, mean)
cdat_av$se <- aggregate(y_norm~SiteID+trait+type+label, cdat2, se)$y_norm
head(cdat_av)

cdat_av$trend <- ifelse(cdat_av$trait %in% c("branch_width","polyp_size", "chl_a","biomass","zoox"), "Increasing","Decreasing")

cdat_av$y_norm.b <- cdat_av$y_norm
cdat_av$y_norm.b[cdat_av$SiteID=="A" & cdat_av$label=="BW"]<-0.19

ggplot(cdat_av, aes(SiteID, y_norm, col=type))+
geom_line(aes(group=trait))+
geom_point(size=0.5)+
geom_path(aes(group=trait))+
geom_errorbar(aes(SiteID, ymin=y_norm-se, ymax=y_norm+se, group=trait), width=0, show.legend = FALSE)+
facet_wrap(~trend, nrow=1)

####################################
####################################
####################################
# Figure 3...
source("figs/fig3.R")
fig3
#ggsave("figs/fig.3b.jpg", fig3, width=145, height=68, units="mm", dpi = 1000)

####################################
####################################
####################################
# figure 4
source("figs/fig4.R")
fig.4
#ggsave("figs/fig.4.jpg", fig.4, width=135, height=110, units="mm", dpi = 1000)

####################################
####################################
####################################
# Effect sizes for fixed and specific CWMs
# 0/NA # lmer/lm # esize/slope

cdat$SiteID2 <- ifelse(cdat$SiteID=="A",1, 
 ifelse(cdat$SiteID=="B",2, 
 ifelse(cdat$SiteID=="C",3, 
 ifelse(cdat$SiteID=="D",4, 
 ifelse(cdat$SiteID=="E",5, 
 ifelse(cdat$SiteID=="F",1, 
 ifelse(cdat$SiteID=="G",2, 
 ifelse(cdat$SiteID=="H",3, 
NA))))))))

params3 <- NULL
table2 <- NULL
clist <- c("specific", "fixed")
for (t in traits){
	for (c in clist){
	df <- subset(cdat, cwm==c)
	df$y <- scale(log(df[,t]))
	df$x <- ifelse(df$Location=="Inner", 1,ifelse(df$Location=="Outer",0,0))
	df <- df[!is.na(df$x),]
	df <- df[!is.na(df$y),]
	mod <-lmer(y ~ x + (1|SiteID2), data=df) 
    coefs <- fixef(mod)
	rsq.m <- r.squaredGLMM(mod)[1]
	rsq.c <- r.squaredGLMM(mod)[2]
	co <- coef(summary(mod))
	se.x <- co["x","Std. Error"]
	pval.x <- co["x","Pr(>|t|)"]
	T.x <- co["x","t value"]
	df.x <- co["x","df"]
	int <- co["(Intercept)","Estimate"]
	se.i <- co["(Intercept)","Std. Error"]
	pval.i <-co["(Intercept)","Pr(>|t|)"]
	T.i <- co["(Intercept)","t value"]
	df.i <- co["(Intercept)","df"]
	#rsq <- r.squaredGLMM(mod)
	se.x <- coef(summary(mod))["x","Std. Error"]
	pval.x <- coef(summary(mod))["x","Pr(>|t|)"]
	params3 <- rbind(params3, data.frame(x=coefs[2], se.x, pval.x, sd.resid=sd(residuals(mod, type="response")), t=t, row.names=NULL, cwm=c))
	table2 <- rbind(table2, data.frame(t=t, cwm=c, E=c("Int", "Slp"), Est=c(int, coefs[2]), se=c(se.i, se.x), df=c(df.i, df.x), T=c(T.i, T.x), P=c(pval.i, pval.x), rsq.m, rsq.c))
	}}

val <- "e.size"
params3$e.size <- params3$x/params3$sd.resid
params3$e.se <- params3$se.x/params3$sd.resid *1.96 
params3$sig <- ifelse(params3$pval.x < 0.05, '*', "")
params3$val <- params3[,val]
params3

pdat <- dcast(t~cwm, data=params3, value.var=val)
pdat$label<- info$label[match(pdat$t, info$trait)]
pdat$type <- info$type[match(pdat$t, info$trait)]
pdat$name <- info$name[match(pdat$t, info$trait)]
pdat$intra <- pdat$specific - pdat$fixed

ggplot()+
geom_bar(data=pdat, aes(fixed, reorder(name, -specific), fill=type), stat="identity", size=0.1, col="black")+
geom_segment(data=pdat, aes(x=fixed, xend=specific, y=reorder(name, -specific),yend=reorder(name, -specific)), stat="identity", fill=NA, size=0.75, col="black", arrow=arrow(length=unit(1,"mm")))+
labs(x="Effect size")+
theme_classic()

pdat$perc <- pdat$intra/pdat$fixed*100
pdat$fac <- pdat$specific/pdat$fixed
pdat

####################################
####################################
####################################

# Sum of sq decomposition withn vs between 
# INTRASPECIFIC = SPECIFIC-FIXED

SSanova <- function(mod){
	SSloc <- summary(mod)[[1]][1,"Sum Sq"]
	SSerr <- summary(mod)[[1]][2,"Sum Sq"]
	SStot <- sum(summary(mod)[[1]][,2])
	Ploc <- summary(mod)[[1]][1,"Pr(>F)"]	
	c(SSloc, SSerr, SStot, Ploc)
}

SSdat <- NULL
Ploc <- NULL
for(t in traits){
#t <- "branch_density"
s.cwm$t <- s.cwm[,t]
f.cwm$t <- f.cwm[,t]
df <- data.frame(ID=rownames(s.cwm), specific=s.cwm$t)
df$fixed <- f.cwm$t[match(df$ID, rownames(f.cwm))]
df$ISV <- df$specific - df$fixed
df$Location <- abun$Location[match(df$ID, abun$ID)]
df$Location <- ifelse(df$Location=="Lagoon", "Outer", as.character(df$Location))
f.aov <- aov(fixed~Location, df) # anovas
s.aov <- aov(specific~Location, df)
i.aov <- aov(ISV~Location, df)
# sum of squares 
temp <-data.frame(trait=t, SS=c("SSloc","SSerr","SStot"), fixed=SSanova(f.aov)[1:3], specific=SSanova(s.aov)[1:3], intra=SSanova(i.aov)[1:3])
temp$covar <- temp$specific - temp$intra - temp$fixed
Pvals <- data.frame(trait=t, fixed=SSanova(f.aov)[4], specific=SSanova(s.aov)[4], intra=SSanova(i.aov)[4])
# proportions
temp$esize <- temp$specific/temp$specific[temp$SS=="SStot"]
#temp$esize <- temp$esize/(1-temp$esize) # cohens f
temp$Prop_turn <- temp$fixed/temp$specific[temp$SS=="SStot"]
temp$Prop_intra <- temp$intra/temp$specific[temp$SS=="SStot"]
temp$Prop_cov <- temp$covar/temp$specific[temp$SS=="SStot"]
temp$Total <- temp$Prop_turn+temp$Prop_intra+temp$Prop_cov
SSdat <- rbind(SSdat, temp)
Ploc <- rbind(Ploc, Pvals)
}

SSdat$type <- info$type[match(SSdat$trait, info$trait)]
SSdat$label <- info$label[match(SSdat$trait, info$trait)]
head(SSdat)

ord <- c("ZD","CH","TB","PB","SD","SV","BD","BW","CW")
SSlong <- melt(SSdat[,c("Prop_intra","Prop_turn","Prop_cov","esize", "type", "label", "SS")])
SSlong$label <- factor(SSlong$label, levels=ord)

ggplot()+geom_bar(data=SSlong, aes(label, value, fill=type), stat="identity")+facet_grid(variable~SS, scales="free")+
geom_hline(yintercept=0)+
labs(x="trait", y="variance explained")

####################################
####################################
####################################
# figure 5

source("figs/fig5.R")
fig.5
#ggsave("figs/fig.5.jpg", fig.5, width=166, height=62, units="mm", dpi = 1000)
