

#######################################
# IMAGES
#######################################
brug<-readPNG("figs/sils/brug.png")
brug<-rasterGrob(brug, interpolate=TRUE)
card<-readPNG("figs/sils/card.png")
card<-rasterGrob(card, interpolate=TRUE)
pdam<-readPNG("figs/sils/pdam.png")
pdam<-rasterGrob(pdam, interpolate=TRUE)
cyl<-readPNG("figs/sils/cyl.png")
cyl<-rasterGrob(cyl, interpolate=TRUE)
gon<-readPNG("figs/sils/gon2.png")
gon<-rasterGrob(gon, interpolate=TRUE)
rus<-readPNG("figs/sils/rus.png")
rus<-rasterGrob(rus, interpolate=TRUE)
stag<-readPNG("figs/sils/stag.png")
stag<-rasterGrob(stag, interpolate=TRUE)

coords <- aggregate(.~Site+SiteID+Location, subset(env, select=-c(Date)), mean)
coords

map2 <- map+ 
   annotation_custom(ggplotGrob(worldmap), xmin=134.17, xmax=134.27, ymin=7.34, ymax=7.398) +
     scalebar(x.min = 134.465, x.max = 134.565, y.min = 7.375, y.max = 7.39, dist = 2, dist_unit = "km",st.size=1.5, st.dist=0.6, box.fill=c("black", "white"),st.bottom = FALSE, transform = TRUE, model = "WGS84", border.size=0.1, height=0.2) +
 geom_point(data=coords, aes(long2, lat), col="black", size=2.3)+
geom_point(data=coords, aes(long2, lat, col=Location), size=2.2)+
scale_colour_manual(values=pal[3:1])+
geom_text(data=coords, aes(long2, lat, label=SiteID), size=1.6, fontface="bold")+
scale_fill_distiller(palette="YlGn", direction=1, breaks=c(29.5, 30), guide=guide_colourbar(ticks.colour = "black", frame.colour="black", title.vjust=0.5,ticks.linewidth = 1, title.position="left", title.angle=90))+
scale_x_continuous(breaks = ewbrks, labels = ewlbls) +
scale_y_continuous(breaks = nsbrks, labels = nslbls) +
guides( col="none")+
  labs(fill="Mean\nSST (°C)")+
  theme_classic()+
  theme(aspect.ratio=0.5, 
  panel.background=element_rect(fill="#78C679"),  
 legend.title=element_text(size=5, face="bold", angle=90), 
 axis.line=element_blank(),
 axis.title=element_blank(),
 legend.background=element_blank(),
 axis.text.y=element_text(size=5, angle=90, hjust=0.5),
 axis.text.x=element_text(size=5), legend.key.height=unit(1.8,"mm"), legend.key.width=unit(1,"mm"), legend.position=c(0.92,0.15),
plot.title=element_text(size=7, face='bold', hjust=0.5),
legend.text=element_text(size=6.5))
map2


####################################
# Environment
envlong$variable2 <- ifelse(envlong$variable=="Temp_deg_C", "Temperature (°C)", 
ifelse(envlong$variable=="Salinity_PSS", "Salinity (S)",
ifelse(envlong$variable=="HDO_mg.l",    "Oxygen (mg/L)",
ifelse(envlong$variable=="pH_units",  "pH" ,  
ifelse(envlong$variable=="Chl_ug.l", "Chlorophyll (μg/L)"  ,  
ifelse(envlong$variable=="Turb_NTU", "Turbidity (NTU)",
as.character(envlong$variable)))))))

envlong$Location <- factor(envlong$Location, levels=c("Outer", "Lagoon", "Inner"))

envplot <-  ggplot(envlong, aes(x=Location, y=value, fill=Location))+
 geom_boxplot(size=0.1, outlier.size=0.05, width=0.5)+
 facet_wrap(~variable2, scales="free_x", nrow=2)+
  scale_fill_manual(values=pal[1:3])+guides(fill="none")+theme_classic()+
 ggtitle("Environmental variation")+
  theme(axis.title.y=element_blank(), 
 panel.spacing.y = unit(1, "mm"),  
 axis.title.x=element_blank(),
 strip.background=element_blank(),
 strip.text=element_text(size=7, margin=margin(0,0,1,0)),
 axis.ticks.length=unit(0.3,"mm"),
 axis.text.x=element_text(size=7, angle=0, hjust=0.5), 
 axis.text.y=element_text(size=7), 
 strip.placement="outside",
 plot.background=element_blank(),
 plot.title=element_text(size=7, face="bold", hjust=0.5),
 panel.background=element_blank(),
 axis.line=element_line(size=0.1))+
 coord_flip()
 envplot

#######################################
# PCA

pcdat <- cbind(na.omit(tdat[,c(traits, "ID","Species", "Genus","SiteID","Location", "Occ")]), data.frame(pc1=pca$x[,1],pc2=pca$x[,2]))
vecs <- data.frame(pc1=pca$rotation[,1],pc2=pca$rotation[,2])

vecs[,colnames(info)] <- info[match(rownames(vecs), info$trait), colnames(info)]
vecs

head(pcdat)
avs <- aggregate(list(pc1 = pcdat$pc1, pc2=pcdat$pc2), by = list(Species = pcdat$Species, Location=pcdat$Location, Occ=pcdat$Occ),  mean)

rectangles2 <- data.frame(type=unique(vecs$type), xmin=0.995, xmax=0.999, ymin=-Inf, ymax=Inf)


traitnames <- ggplot()+geom_text(data=vecs, aes(x=1, y=label, label=paste(name,"-", label)), size=2.3, hjust=0)+
lims(x=c(0.99,1.05))+
geom_rect(data=rectangles2, aes(xmin=xmin, ymin= ymin,ymax= ymax, xmax=xmax, fill=type))+
geom_text(data=rectangles2, aes(x=1, y=c(5,2,5), label=type), size=2.3, hjust=0, fontface="bold")+
facet_grid(rows=vars(type), scales="free", space="free")+
theme_void()+
scale_fill_manual(values=pal2[3:1])+
guides(fill="none")+
coord_cartesian(clip="off")+
theme(panel.spacing=unit(2,"mm"), 
plot.margin = unit(c(1,5,9,-3), "mm"),
strip.text=element_blank())
traitnames

ex <- 3
ex2 <- 4

exp<-round(c(summary(pca)[[1]][1]^2/sum(summary(pca)[[1]]^2),summary(pca)[[1]][2]^2/sum(summary(pca)[[1]]^2)),3)*100
exp

###################

# original vectors

pcaX <- prcomp(na.omit(log(tdat[tdat$Occ %in% c("Both","dunno"),traits])), center=T, scale=T)
biplot(pcaX)
vecs <- data.frame(pc1=pcaX$rotation[,1],pc2=pcaX$rotation[,2])

vecs[,colnames(info)] <- info[match(rownames(vecs), info$trait), colnames(info)]
vecs

vecs$pc1b <-ifelse(vecs$label=="ZD", vecs$pc1-0.03,vecs$pc1)
vecs$pc2b <-ifelse(vecs$label=="CH", vecs$pc2+0.03,vecs$pc2)
vecs$pc1b[vecs$label=="CH"] <-vecs$pc1b[vecs$label=="CH"]+0.03

pcplot1 <- ggplot(pcdat, aes(pc1, pc2))+
geom_segment(data=vecs, aes(pc1*ex, pc2*ex, yend=0, xend=0), size=0.3, col="black")+
geom_point(data=vecs, aes(pc1b*ex2, pc2b*ex2), size=4.5, col="white")+
geom_point(data=vecs, aes(pc1b*ex2, pc2b*ex2, col=type), size=4.5, alpha=0.5)+
geom_point(data=vecs, aes(pc1b*ex2, pc2b*ex2, col=type), size=4.5, shape=21, stroke=0.45)+
geom_text(data=vecs, aes(pc1b*ex2, pc2b*ex2, label=label),size=2,fontface="bold")+
ggtitle("")+
lims(x=c(min(vecs$pc1b*ex2)*1.1, max(vecs$pc1*ex2)*1.1), y=c(min(vecs$pc2b*ex2)*1.1, max(vecs$pc2b*ex2)*1.1))+
scale_colour_manual(values=pal2[3:1])+
guides(fill="none", col="none")+
theme_bw()+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), 
axis.title=element_text(size=8),
legend.position=c(0.8,0.2),
legend.background=element_blank(),
axis.text=element_text(size=8),
plot.title=element_text(size=9, hjust=0.5, face="bold"),
legend.title=element_blank(), legend.text=element_text(size=7), legend.key.size=unit(0.1,"mm"), legend.direction="vertical")
pcplot1 

a_loc <- aggregate(cbind(n=a_tran$n, cover=a_tran$cover)~Species+Location, a_tran, mean)

# calculate relative abundance
tot <- aggregate(cover~Transect+SiteID+Circle, a_tran, sum)

avs$spp <- spp$spp[match(avs$Species, spp$species)]
avsmax <- aggregate(pc2~spp, avs, max)
avsmax$pc1 <- avs$pc1[match(avsmax$pc2, avs$pc2)]
pcdat$cover <- a_loc$cover[match(paste(pcdat$Species, pcdat$Location), paste(a_loc$Species, a_loc$Location))]
tot <- aggregate(cover~Location, a_loc, sum)
pcdat$tot <- tot$cover[match(pcdat$Location, tot$Location)]
pcdat$relabun <- pcdat$cover/pcdat$tot
avsmax


avsmax$pc1b <- avsmax$pc1+0.2
avsmax$pc2b <- avsmax$pc2+0.2

avsmax$pc1b[2] <- avsmax$pc1b[2]-1
avsmax$pc2b[2] <- avsmax$pc2b[2]-0.7

avsmax$pc1b[1] <- avsmax$pc1b[1]+0.4
avsmax$pc2b[1] <- avsmax$pc2b[1]+0.1

avsmax$pc1b[3] <- avsmax$pc1b[3]+0.3
avsmax$pc2b[3] <- avsmax$pc2b[3]-0.4

avsmax$pc1b[4] <- avsmax$pc1b[4]-0.7

avsmax$pc1b[5] <- avsmax$pc1b[5]+0.5
avsmax$pc2b[5] <- avsmax$pc2b[5]-0.4

avsmax$pc1b[6] <- avsmax$pc1b[6]-0.5
avsmax$pc2b[6] <- avsmax$pc2b[6]-0.05

avsmax$pc1b[7] <- avsmax$pc1b[7]+0.2
avsmax$pc2b[7] <- avsmax$pc2b[7]-0.2


pcplot2 <- ggplot(pcdat[pcdat$Occ %in% c("Both", "dunno"),], aes(pc1, pc2))+
annotation_custom(gon, xmin=2.1, xmax=4, ymin=2.1, ymax=3)+
annotation_custom(stag, xmin=-0.4, xmax=1, ymin=2.5, ymax=3.5)+
annotation_custom(brug, xmin=0.2, xmax=1.6, ymin=1.1, ymax=3.1)+
annotation_custom(pdam, xmin=-1.3, xmax=0, ymin=2.8, ymax=3.7)+
annotation_custom(card, xmin=-2.5, xmax=-1, ymin=2.7, ymax=3.8)+
annotation_custom(rus, xmin=1.2, xmax=2.8, ymin=0.5, ymax=2.2)+
annotation_custom(cyl, xmin=1.5, xmax=2.5, ymin=-0.7, ymax=1)+
#annotation_custom(ggplotGrob(occboth), xmin=3.5, xmax=6, ymin=-4.4, ymax=-1.6)+
geom_point(aes(col=Location, size=relabun), alpha=0.75)+
#scale_colour_manual(values=colsT2)+
scale_colour_manual(values=pal[1:3])+
labs(x=paste("PC1 (", exp[1],"%)", sep=""), y=paste("PC2 (", exp[2],"%)", sep=""))+
geom_path(data=avs[avs$Occ %in% c("Both", "dunno"),], aes(pc1, pc2, group=Species), arrow=arrow(length=unit(1, "mm"),  type="closed"), size=0.4)+
#geom_label_repel(data=avsmax, aes(pc1+0.2,pc2+0.2, label=spp),  col="slategrey", force=0.001, fontface="bold", label.padding=0, label.size=0, alpha=0.7, size=2)+
geom_text(data=avsmax, aes(pc1b,pc2b, label=spp),  col="slategrey", fontface="bold",  size=2.3)+
xlim(c(-4,5.5))+ylim(c(-4,4))+
ggtitle("Intraspecific trait change")+
#guides(fill="none", col="none")+
scale_radius(range=c(1,2.5))+
guides(size="none")+
theme_bw()+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), 
axis.title=element_text(size=7),
axis.text=element_text(size=7),
plot.title=element_text(size=8, hjust=0.5, face="bold"),
legend.title=element_blank(), 
legend.position="top",
plot.background=element_blank(),
legend.spacing=unit(0.5,"mm"),
legend.direction="horizontal",
legend.margin=margin(0, 0, -10, 0),
legend.text=element_text(size=7), legend.key.size=unit(0.2,"mm"))
pcplot2


fig.1 <- plot_grid(
plot_grid(
plot_grid(NULL,map2+ggtitle("Sampling sites"), rel_widths=c(0.03,1)),
envplot+facet_wrap(~variable2, nrow=3, scales="free_x")
, nrow=2, rel_heights=c(1,1.3), labels=c("A","B"), label_size=9, vjust=3),
plot_grid(pcplot2,  plot_grid(NULL,
plot_grid(NULL,pcplot1+theme_void(),NULL, ncol=1, rel_heights=c(-0.2,1,0.2)),
plot_grid(traitnames, NULL, ncol=1,rel_heights=c(1, -0.1)),
rel_widths=c(0.2,1,1), nrow=1),
ncol=1, rel_heights=c(1,0.5), labels=c("C","D"), label_size=9, vjust=c(3,1))
, rel_widths=c(1.1,1))
fig.1

