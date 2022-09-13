

# Supplementary figure 1

#library("psych")
#pairs.panels(na.omit(log(tdat[,traits])), scale=T, cex.cor=2)


####################################
####################################
####################################

# Supplementary figure 2


tlong <- melt(dat[,c("Species","Location","Genus","SiteID","Occ", traits)], id.var=c("SiteID","Location","Species", "Genus", "Occ"), variable.name="trait")
head(tlong)

tlong$Location2 <- ifelse(tlong$Location=="Lagoon", "Outer", as.character(tlong$Location))

head(tlong)
spmean<-aggregate(value~trait+Species, tlong, mean)
tlong$spmean <- spmean$value[match(paste(tlong$t, tlong$Species), paste(spmean$t, spmean$Species))]
tlong$value2 <- tlong$value/tlong$spmean
tlong$name <- info$name[match(tlong$trait, info$trait)]
head(tlong)


ggplot(tlong[tlong$Occ=="Both",], aes(SiteID, value))+
geom_rect(data=NULL, aes(xmin="F", xmax="H", ymin=-Inf, ymax=Inf), fill="grey90")+
stat_summary(fun = mean, geom = "path", aes(col=Species, group=Species), size=0.25)+
stat_summary(fun = mean, geom = "point", aes(col=Species), size=0.7)+
stat_summary(fun.data = mean_se, geom = "errorbar", aes(col=Species), width=0, size=0.25)+
facet_wrap(~name, scales="free_y", ncol=3)+
#scale_y_log10()+
scale_colour_manual(values=colsT2[names(colsT2) %in% unique(tlong[tlong$Occ=="Both","Species"])])+
labs(y="Trait value", x="Site")+
theme_classic()+theme(strip.background=element_blank(), axis.line=element_line(size=0.1), legend.text=element_text(face="italic", size=8), legend.position="bottom", legend.title=element_blank(), legend.key.width=unit(1,"mm"))


#######################

# Supplementary figure 3

head(s.cwm)
head(cdat)
head(dat)

head(cdat)
udat <- cdat[cdat$cwm=="specific",]
udat2 <- cdat[cdat$cwm=="fixed",]
plots <- list()
for (t in traits){
udat$t <- udat[,t]
udat2$t <- udat2[,t]
dat$t <- dat[,t]
plots[[t]] <- ggplot()+
#geom_point(data=udat, aes(SiteID, t))+
stat_summary(data=dat,fun="mean", geom="line",aes(SiteID, t, group=spp), col="grey")+
ggtitle(t)+
stat_summary(data=udat,fun="mean", geom="line",aes(SiteID, t, group=1))+
stat_summary(data=udat,fun="mean", geom="line",aes(SiteID, t, group=1))+
stat_summary(data=udat2,fun="mean", geom="line",aes(SiteID, t, group=1), col="darkred")+
theme_classic()}
plot_grid(plotlist=plots)




#######################

# Supplementary figure 4

SSlong$variableX <- ifelse(SSlong$variable=="esize","Effect size",ifelse(SSlong$variable=="Prop_intra","Intraspecific", ifelse(SSlong$variable =="Prop_turn","Species composition",ifelse(SSlong$variable=="Prop_cov","Covariation", NA))))
SSlong$variableX <- factor(SSlong$variableX, levels=c("Effect size","Intraspecific","Species composition","Covariation"))

ggplot()+geom_bar(data=SSlong, aes(label, value, fill=type), stat="identity")+facet_grid(variableX~SS, scales="free")+
scale_fill_manual(values=pal2[3:1])+
theme_classic()+theme(strip.text=element_text(), strip.background=element_blank(), legend.title=element_blank())+
geom_hline(yintercept=0)+
labs(x="trait", y="variance explained")

