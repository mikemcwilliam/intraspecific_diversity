

####################################
# Sampling/composition
lines <- data.frame(SiteID = unique(env$SiteID))
lines$Location <- env$Location[match(lines$SiteID, env$SiteID)]
lines$down<-ifelse(lines$SiteID %in% c("B","G"),"n","y")
lines$text<-ifelse(lines$SiteID %in% c("B","G","D"),"y","n")

lineplot <- ggplot()+
geom_line(data=lines, aes(x=SiteID, y=1, col=Location, group=Location))+
geom_text(data=lines[lines$text=="y",], aes(x=SiteID, y=2, label=Location), size=2.5, nudge_x=c(0,0.5,0))+
lims(y=c(0,3))+
scale_x_discrete(expand=c(0.01,0.01))+
geom_segment(data=lines[lines$down=="y",], aes(x=SiteID, xend=SiteID, y=1, yend=0.5, col=Location))+
guides(col="none", fill="none")+
ggtitle("Species abundances")+
scale_colour_manual(values=pal[3:1])+
theme_void()+
theme(plot.title=element_text(size=8, hjust=0.5, face="bold"))
lineplot

orderT2 <- c("Porites cf. cylindrica" ,"Porites cf. rus" , "Porites cf. nigrescens",  "Porites Massive",
"Pocillopora cf. damicornis", "Pocillopora cf. verrucosa",    "Stylophora cf. pistillata" ,
"Goniastrea spp." ,   "Lobophyllia spp.",
"Anacropora spinosa", 
"Isopora cf. brueggemanni",  
"Acropora cf. tenuis","Acropora cf. muricata","Acropora cf. subglabra", "Acropora cf. hyacinthus", "Acropora cf. humilis")

colsT2 <- as.character(c("#762a83","#9970ab","#c2a5cf","#e7d4e8", 
"#d95f02","#fee0b6","#e6ab02",
"#5aae61","#1b7837",
"#f1b6da",
"#c51b7d",
"#313695","#4575b4","#abd9e9","#74add1", "#23545c"))
names(colsT2) <- orderT2

unique(a_site$Species)

a_site$Species <- factor(a_site$Species, levels=c(orderT2, UNSAMP))
a_site$SpeciesS <- factor(a_site$SpeciesS, levels=c(orderT2, "Other"))

a_site$spp <- as.character(spp$spp[match(a_site$SpeciesS,spp$species)])

#a_site$spp[a_site$spp=="Gon"]<-"Gre/Gfa"
a_site$spp[a_site$spp=="Lob"]<-"Lhe/Lco"

a_site$SpeciesS2 <- paste(a_site$SpeciesS," (", a_site$spp,")", sep="")
a_site$SpeciesS2[a_site$SpeciesS2=="Other (NA)"]<-"Other"

colsT3 <- colsT2
names(colsT3) <- a_site$SpeciesS2[match(names(colsT2),a_site$SpeciesS)]
colsT3
colsT3 <- c(colsT3, Other="Grey")

a_site$SpeciesS2 <- factor(a_site$SpeciesS2, levels=names(colsT3))


strato <- ggplot(a_site, aes(x=SiteID, y=cover, group=Species, fill=SpeciesS2))+ 
geom_area(position = 'stack', col="black", size=0.1)+
scale_x_discrete(expand=c(0,0))+
scale_y_continuous(expand=c(0,0))+
labs(y="% cover", x="Site")+
theme_classic()+
guides(fill=guide_legend(ncol=2))+
scale_fill_manual(values=colsT3)+
#ggtitle("Species abundances")+
theme(legend.title=element_blank(), 
#axis.title.x=element_blank(),
legend.position="bottom",
axis.title=element_text(size=8),
axis.text=element_text(size=8),
legend.text=element_text(size=7, face="italic"), 
plot.title=element_text(size=8, face="bold", hjust=0.5),
legend.key.size=unit(2,"mm"))
strato

#plot_grid(lineplot, strato,align="v", axis="lr", rel_heights=c(0.07, 1), ncol=1)


legendSp <- get_legend(strato+theme(legend.background=element_blank()))

samp_plot <- plot_grid(lineplot, strato+guides(fill="none"),  ncol=1, rel_heights=c(0.2, 1), align="v", axis="lr")
samp_plot 

####################################
####################################
####################################
# Locally dominant species 

pcdat <- cbind(na.omit(dat[,c(traits, "ID","Species","Species_TB", "Genus","SiteID","Location", "Occ", "spp")]), data.frame(pc1=pca3$x[,1],pc2=pca3$x[,2]))

vecs <- data.frame(pc1=pca3$rotation[,1],pc2=pca3$rotation[,2])
vecs[,colnames(info)] <- info[match(rownames(vecs), info$trait), colnames(info)]
vecs

ggplot(pcdat, aes(pc1, pc2))+geom_point(aes(col=Occ))



vecs$pc1b <-ifelse(vecs$label=="ZD", vecs$pc1-0.03,vecs$pc1)
vecs$pc2b <-ifelse(vecs$label=="CH", vecs$pc2+0.03,vecs$pc2)
vecs$pc1b[vecs$label=="CH"] <-vecs$pc1b[vecs$label=="CH"]+0.03

ex3 <- 2
ex4 <- 3.5

pcplotV <- ggplot(pcdat, aes(pc1, pc2))+
geom_segment(data=vecs, aes(pc1*ex3, pc2*ex3, yend=0, xend=0), size=0.2, col="black")+
geom_point(data=vecs, aes(pc1b*ex3, pc2b*ex3), size=1.7, col="white")+
geom_point(data=vecs, aes(pc1b*ex3, pc2b*ex3, col=type), size=1.7, alpha=0.8)+
geom_point(data=vecs, aes(pc1b*ex3, pc2b*ex3), size=1.7, shape=21, stroke=0.1)+
#geom_text_repel(data=vecs, aes(pc1b*ex4, pc2b*ex4, label=label),size=1.5,fontface="bold", force=0.006)+
lims(x=c(min(vecs$pc1b*ex3)*1.2, max(vecs$pc1*ex3)*1.2), y=c(min(vecs$pc2b*ex3)*1.15, max(vecs$pc2b*ex3)*1.15))+
scale_colour_manual(values=pal2[3:1])+
guides(fill="none", col="none")+
theme_void()+theme(panel.border=element_rect(fill=NA, colour="white"))
pcplotV

pcplotV2 <- ggplot(pcdat, aes(pc1, pc2))+
geom_segment(data=vecs, aes(xend=pc1*ex3, yend=pc2*ex3, y=0, x=0, col=type), size=0.3, arrow=arrow(type="closed", length=unit(0.5,"mm")))+
lims(x=c(min(vecs$pc1b*ex3)*1.2, max(vecs$pc1*ex3)*1.2), y=c(min(vecs$pc2b*ex3)*1.15, max(vecs$pc2b*ex3)*1.15))+
scale_colour_manual(values=pal2[3:1])+
guides(fill="none", col="none")+
theme_void()+theme(panel.border=element_rect(fill=NA, colour="grey"))
pcplotV2


sp.avs <- aggregate(pc1~spp+Location+Occ, pcdat, mean)
sp.avs$pc2 <- aggregate(pc2~spp+Location+Occ, pcdat, mean)$pc2
sp.avs <- sp.avs[!sp.avs$Occ =="Both",]
sp.avs

sp.avs2 <- sp.avs[!sp.avs$Location=="Lagoon",]
sp.avs2$pc1b <- sp.avs2$pc1 - 1.3
sp.avs2$pc2b <- sp.avs2$pc2 + 0
sp.avs2$pc1b[4] <- sp.avs2$pc1b[4] + 2.6
sp.avs2$pc1b[3] <- sp.avs2$pc1b[3] + 2.6
sp.avs2$pc2b[3] <- sp.avs2$pc2b[3] -0.2
sp.avs2$pc1b[8] <- sp.avs2$pc1b[8] + 2.6

exp<-round(c(summary(pca3)[[1]][1]^2/sum(summary(pca3)[[1]]^2),summary(pca3)[[1]][2]^2/sum(summary(pca3)[[1]]^2)),3)*100
exp

xlab1 <- ggplot(data=NULL)+
geom_text(aes(0.9,1,label="Branching"), col="grey", size=2.5)+
geom_text(aes(2.1,1,label="Massive"), col="grey", size=2.5)+
geom_segment(aes(0.35,1,xend=0, yend=1), arrow=arrow(length=unit(1,"mm")), col="grey")+
geom_segment(aes(2.65,1,xend=3, yend=1), arrow=arrow(length=unit(1,"mm")), col="grey")+
theme_void()
xlab1

ylab1 <- ggplot(data=NULL)+
geom_text(aes(1,1,label="Tissue"), col="grey", size=2.5, angle=270)+
geom_text(aes(2,1,label="Skeleton"), col="grey", size=2.5, angle=270)+
geom_segment(aes(0.5,1,xend=0, yend=1), arrow=arrow(length=unit(1,"mm")), col="grey")+
geom_segment(aes(2.5,1,xend=3, yend=1), arrow=arrow(length=unit(1,"mm")), col="grey")+
theme_void()+coord_flip()
ylab1

domplot <- ggplot(pcdat, aes(pc1, pc2))+
#annotation_custom(anac, xmin=-5, xmax=-3, ymin=-3.2, ymax=-2.2)+
#annotation_custom(lobo, xmin=4, xmax=8, ymin=-3.3, ymax=-2.3)+
annotation_custom(ggplotGrob(pcplotV2), xmin=4.5, xmax=7.8, ymin=1.5, ymax=3.2)+
#annotation_custom(ggplotGrob(xlab1), xmin=-5, xmax=7.5, ymin=3.2, ymax=3.5)+
#annotation_custom(ggplotGrob(ylab1), xmin=7, xmax=7.7, ymin=-2.5, ymax=2.5)+
#geom_point(data=pcdat[!pcdat$Occ %in% c("Both","dunno"),],aes(col=Location), size=2)+
geom_point(col="grey85", size=0.1)+
#geom_point(data=pcdat[!pcdat$Occ=="Both",],aes(col=Location), size=0.5, shape=21, stroke=0.25)+
#geom_segment(data=vecs, aes(0,0, xend=pc1*4.5, yend=pc2*4.5), size=0.1)+
#geom_text_repel(data=vecs,aes(pc1*5.5, pc2*5.5, label=label), size=2.5, col="slategrey", force=0.0025)+
geom_point(data=sp.avs, aes(col=Location), size=3, alpha=0.75)+
geom_path(data=sp.avs, aes(group=spp), size=0.2)+
geom_text(data=sp.avs2, aes(x=pc1b,y=pc2b, label=spp), size=2, fontface="bold", col="slategrey")+
labs(x=paste("PC1 (", exp[1],"%)", sep=""), y=paste("PC2 (", exp[2],"%)", sep=""))+
scale_colour_manual(values=pal[1:3])+
#ggtitle("Species turnover")+
guides(col="none",fill="none")+
lims(y=c(-3.3,3))+
theme_bw()+
scale_x_continuous(breaks=c(-3,0,3,6), limits=c(-5,7.5))+
theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
axis.text=element_text(size=8),
axis.title=element_text(size=8),
plot.title=element_text(size=8, face="bold", hjust=0.5))
domplot 


pcdatX <- pcdat
pcdatX$Species_TB[pcdatX$Species_TB=="Porites attenuata"]<-"Porites cf. cylindrica"
pcdatX$sppX <- spp$spp[match(pcdatX$Species_TB, spp$species)]
#pcdatX$Occ[pcdatX$sppX=="Gfa"]<-"Inner"
#pcdatX$Occ[pcdatX$sppX=="Gre"]<-"Outer"
#pcdatX$Occ[pcdatX$sppX=="Gst"]<-"Inner"

pcdatX<-pcdatX[!pcdatX$ID=="s2_goni_1",]

sp.avsX <- aggregate(pc1~sppX+Location+Occ, pcdatX, mean)
sp.avsX$pc2 <- aggregate(pc2~sppX+Location+Occ, pcdatX, mean)$pc2
sp.avsX <- sp.avsX[!sp.avsX$Occ =="Both",]
sp.avsX

#sp.avs2X <- sp.avsX[!sp.avsX$Location=="Lagoon",]
#sp.avs2X
#sp.avs2X$pc1b <- sp.avs2X$pc1 - 1.3
#sp.avs2X$pc2b <- sp.avs2X$pc2 + 0
#sp.avs2X$pc1b[1] <- sp.avs2X$pc1b[1] + 2.6
#sp.avs2X$pc1b[5] <- sp.avs2X$pc1b[5] + 1
#sp.avs2X$pc2b[5] <- sp.avs2X$pc2b[5] + 0.5
#sp.avs2X$pc1b[3] <- sp.avs2X$pc1b[3] + 2
#sp.avs2X$pc2b[3] <- sp.avs2X$pc2b[3] -0.4
#sp.avs2X$pc1b[10] <- sp.avs2X$pc1b[10] + 2.6

sp.avs2X <- sp.avsX[!sp.avsX$Location=="Lagoon",]
sp.avs2X$pc1b <- sp.avs2X$pc1 - 1.3
sp.avs2X$pc2b <- sp.avs2X$pc2 + 0
sp.avs2X$pc1b[2] <- sp.avs2X$pc1b[2] + 2.6
sp.avs2X$pc1b[3] <- sp.avs2X$pc1b[3] + 2.6
sp.avs2X$pc1b[4] <- sp.avs2X$pc1b[4] - 0.5
sp.avs2X$pc2b[5] <- sp.avs2X$pc2b[5] + 0.3

sp.avs2X$pc1b[6] <- sp.avs2X$pc1b[6] + 1.2
sp.avs2X$pc2b[6] <- sp.avs2X$pc2b[6] - 0.4

sp.avs2X$pc1b[7] <- sp.avs2X$pc1b[7] + 1
sp.avs2X$pc2b[7] <- sp.avs2X$pc2b[7] + 0.4
######################
# CWM of pca axes??
pc.cwmX<-aggregate(pc1~Location, cdat[cdat$cwm=="specific",], mean)
pc.cwmX$pc2<-aggregate(pc2~Location, cdat[cdat$cwm=="specific",], mean)$pc2
#pc.cwmX<-aggregate(pc1~Location, cdat[cdat$cwm=="fixed",], mean)
#pc.cwmX$pc2<-aggregate(pc2~Location, cdat[cdat$cwm=="fixed",], mean)$pc2
pc.cwmX


domplotX <- ggplot(pcdatX, aes(pc1, pc2))+
#annotation_custom(anac, xmin=-5, xmax=-3, ymin=-3.2, ymax=-2.2)+
#annotation_custom(lobo, xmin=4, xmax=8, ymin=-3.3, ymax=-2.3)+
annotation_custom(ggplotGrob(pcplotV2), xmin=4.5, xmax=7.8, ymin=1.5, ymax=3.2)+
#annotation_custom(ggplotGrob(xlab1), xmin=-5, xmax=7.5, ymin=3.2, ymax=3.5)+
#annotation_custom(ggplotGrob(ylab1), xmin=7, xmax=7.7, ymin=-2.5, ymax=2.5)+
#geom_point(aes(col=Location), size=2)+
geom_point(col="grey85", size=0.1)+
#geom_point(data=pcdatX[!pcdatX$Occ %in% c("Both","dunno"),],aes(col=Location), size=2)+
#geom_point(data=pcdat[!pcdat$Occ=="Both",],aes(col=Location), size=0.5, shape=21, stroke=0.25)+
#geom_segment(data=vecs, aes(0,0, xend=pc1*4.5, yend=pc2*4.5), size=0.1)+
geom_path(data=pc.cwmX, aes(pc1, pc2),  arrow=arrow(length=unit(1, "mm"),  type="closed", end="first"), size=0.4)+
#geom_text_repel(data=vecs,aes(pc1*5.5, pc2*5.5, label=label), size=2.5, col="slategrey", force=0.0025)+
geom_point(data=sp.avsX, aes(col=Location), size=3, alpha=0.75)+
geom_path(data=sp.avsX, aes(group=sppX), size=0.2, col="grey50")+
geom_text(data=sp.avs2X, aes(x=pc1b,y=pc2b, label=sppX), size=2, fontface="bold", col="slategrey")+
labs(x=paste("PC1 (", exp[1],"%)", sep=""), y=paste("PC2 (", exp[2],"%)", sep=""))+
scale_colour_manual(values=pal[1:3])+
#ggtitle("Species turnover")+
guides(col="none",fill="none")+
lims(y=c(-3.3,3))+
theme_bw()+
scale_x_continuous(breaks=c(-3,0,3,6), limits=c(-5,7.5))+
theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
axis.text=element_text(size=8),
axis.title=element_text(size=8),
plot.title=element_text(size=8, face="bold", hjust=0.5))
domplotX





domplot2 <- plot_grid(plot_grid(NULL, domplotX, nrow=2, rel_heights=c(0.1,1)), NULL, nrow=1, rel_widths=c(1,0.07))+
draw_grob(ggplotGrob(ylab1), 0.88, 0.2, width=0.14, height=0.65)+
draw_grob(ggplotGrob(xlab1), 0.2, 0.85, width=0.65, height=0.13)+
draw_label("Species turnover", 0.5, 0.97, size=8, fontface="bold")
domplot2


cwmplot<-ggplot(cdat_av, aes(SiteID, y_norm, col=type))+
geom_line(aes(group=trait))+
geom_point(size=0.5)+
geom_path(aes(group=trait))+
geom_errorbar(aes(SiteID, ymin=y_norm-se, ymax=y_norm+se, group=trait), width=0, show.legend = FALSE)+
geom_label(data=cdat_av[cdat_av$SiteID=="H" & cdat_av$trend=="Decreasing",], aes(label=label, fill=type), col="black", nudge_x=1, size=2, alpha=0.4, label.padding=unit(0.5,"mm"), label.size=0, fontface="bold", show.legend = FALSE)+
geom_label(data=cdat_av[cdat_av$SiteID=="A" & cdat_av$trend=="Increasing",], aes(y=y_norm.b, label=label, fill=type), col="black", nudge_x=-1, size=2, alpha=0.4, label.padding=unit(0.5,"mm"), label.size=0, fontface="bold", show.legend = FALSE)+
facet_wrap(~trend, nrow=1)+
guides(col = guide_legend(direction="horizontal", override.aes = list(size = 2)))+
labs(x="Site", y="Community trait mean (scaled)")+
#scale_colour_manual(values=c("red", "black"))+
scale_colour_manual(values=pal2[3:1])+
scale_fill_manual(values=pal2[3:1])+
theme_classic()+
ggtitle("Community trait change")+
expand_limits(x= c(-0.75, 10))+
scale_y_continuous(breaks=c(0,0.5,1))+
theme(strip.background=element_blank(), 
strip.text=element_text(size=8),
legend.title=element_blank(),
legend.text=element_text(size=8),
legend.margin=margin(-1,0,-7,0),
legend.key.size=unit(1,"mm"),
legend.position="top",
axis.line=element_line(size=0.1), 
plot.title=element_text(size=8, face="bold", hjust=0.5),
axis.text=element_text(size=8), axis.title=element_text(size=8))
cwmplot


fig.4 <- plot_grid(
plot_grid(domplot2, 
plot_grid(samp_plot,NULL,rel_widths=c(1,0.05)), ncol=1, labels=c("A", "B"), label_size=9,rel_heights=c(1,1)), 
plot_grid(cwmplot, NULL,legendSp,NULL,
 ncol=1, rel_heights=c(1,-0.05, 0.41,0.05), labels=c("C","",""), label_size=9), rel_widths=c(1,1.5))+
draw_line(x=c(0.37,0.41), y=c(0.26,0.26), linetype="dotted")
fig.4




