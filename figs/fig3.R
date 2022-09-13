
# factors 
mus <- aggregate(zoox~SiteID+Species, tdat[tdat$Occ=="Both",], mean)
factors<-data.frame(sp=unique(mus$Species), zoox=aggregate(zoox~Species, mus, max)$zoox/aggregate(zoox~Species, mus, min)$zoox)
mus <- aggregate(chl_a~SiteID+Species, tdat[tdat$Occ=="Both",], mean)
factors$chl_a <- aggregate(chl_a~Species, mus, max)$chl_a/aggregate(chl_a~Species, mus, min)$chl_a
factors

ggplot(factors, aes(zoox, chl_a))+geom_text(aes(label=sp))+geom_smooth(method="lm")

means <- aggregate(chl_a~SiteID+Species, tdat[tdat$Occ=="Both",], mean)
Hmeans <- subset(means, SiteID=="H")
Hmeans <- rbind(Hmeans, means[means$Species=="Acropora cf. muricata" & means$SiteID=="F",])
Hmeans <- rbind(Hmeans, means[means$Species=="Isopora cf. brueggemanni" & means$SiteID=="G",])
Hmeans <- rbind(Hmeans, means[means$Species=="Pocillopora cf. damicornis" & means$SiteID=="G",])
Hmeans$spp <- spp$spp[match(Hmeans$Species, spp$species)]
#Hmeans$end <- c(8,8,8,8,8,8,8)
Hmeans$end <- c(8,8,8,8,7,7.5,7.5)

mod <- lm(log(chl_a)~log(zoox), dat)
summary(mod)
confint(mod)

lines <- data.frame(SiteID = unique(env$SiteID))
lines$Location <- env$Location[match(lines$SiteID, env$SiteID)]
lines$down<-ifelse(lines$SiteID %in% c("B","G"),"n","y")
lines$text<-ifelse(lines$SiteID %in% c("B","G","D"),"y","n")

tdat$mID <- paste(tdat$Species, tdat$SiteID)
means$mID <- paste(means$Species, means$SiteID)
tdat$mean <- means$chl_a[match(tdat$mID, means$mID)]

avs <- aggregate(chl_a~Species+SiteID+Location+Occ2, dat, mean)
avs$cov <- a_site$cover[match(paste(avs$SiteID, avs$spp), paste(a_site$SiteID, a_site$spp))]

udat <- cdat[cdat$cwm=="specific",]
lavs <- aggregate(chl_a~Location, udat, mean)
lavs$se <- aggregate(chl_a~Location, udat, se)$chl_a
lavs$min <- aggregate(chl_a~Location, udat, min)$chl_a
lavs$max <- aggregate(chl_a~Location, udat, max)$chl_a

chloroplot1 <- ggplot(tdat[tdat$Occ%in%"Both",], aes(SiteID, chl_a))+
geom_point(data=tdat, aes(SiteID, chl_a), col="white")+
geom_rect(data=lavs[lavs$Location=="Inner",], inherit.aes=F, aes(xmin=9.5, xmax=10, ymin=min, ymax=max), fill=pal[3], alpha=0.5)+
geom_rect(data=lavs[lavs$Location=="Inner",], inherit.aes=F, aes(xmin=9.5, xmax=10, ymin=chl_a-se*2, ymax=chl_a+se*2), fill=pal[3], alpha=1, fill="white")+
geom_segment(data=Hmeans, aes(x=SiteID, xend=end, yend=chl_a, y=chl_a), size=0.2, col="grey")+
geom_line(data=lines, aes(x=SiteID, y=0, col=Location, group=Location), size=2)+
geom_line(aes(SiteID, mean, group=Species),size=0.2, na.rm=T, col=pal2[1])+
geom_point(aes(SiteID, mean, group=Species), shape=21, fill=pal2[1], stroke=0.1,size=1)+
scale_colour_manual(values=pal[3:1])+
scale_y_log10()+
#expand_limits(x= c(0, 9.5))+
guides(col="none", fill="none")+
labs(y=expression(~Chlorophyll~content~"("~µg~cm^-2~")"), x="Site")+
theme_classic()+theme(axis.text=element_text(size=8), axis.title.y=element_text(size=8), 
axis.title.x=element_text(size=8,margin=margin(-5,0,0,0)), 
axis.line=element_line(size=0.1), plot.backgroun=element_blank())
chloroplot1


chloroplot2 <- ggplot(tdat[tdat$Occ %in% c("Both"),], aes(zoox, chl_a))+
 geom_point(aes(col=Location), size=1)+
geom_smooth(aes(col=Location), method="lm",  size=0.75, se=F)+
scale_y_log10()+
scale_x_log10()+
scale_colour_manual(values=pal[1:3])+
scale_shape_manual(values=c(4,16))+
theme_classic()+
ggtitle("")+
labs(x=expression(~Symbiont~density~"("*10^6~cells~cm^-2*")"))+
theme(axis.line=element_line(size=0.1),
axis.text.y=element_text(size=8),
axis.text.x=element_text(size=8),
axis.title.y=element_blank(), 
legend.title=element_blank(),
legend.position=c(0.2, 0.9),
legend.background=element_blank(),
legend.key.size=unit(1, "mm"),
legend.text=element_text(size=8),
plot.title=element_text(size=8, face="bold", hjust=0.5),
axis.title.x=element_text(size=8, margin=margin(-5,0,0,0)))
chloroplot2

chloroplot3 <- ggplot(tdat[tdat$Occ %in% "Both",], aes(x=Location, y=chl_a))+
geom_boxplot(fill="white", aes(col=Location), outlier.shape=NA, width=0.2,position=position_nudge(0.4), coef=0)+
geom_bar(aes(y=Inf), stat="identity", width=0.75, fill="white")+
geom_line(col="black", size=0.2)+
geom_jitter(shape=21, stroke=0.1, aes(fill=ID), height=0, width=0.15, size=3)+
guides(fill="none", col="none")+
scale_colour_manual(values=pal[1:3])+
scale_fill_manual(values=shades)+
scale_y_log10()+
theme_classic()+
scale_x_discrete(expand=c(0,1))+
theme(axis.line=element_line(size=0.1),
axis.text.y=element_text(size=8),
axis.text.x=element_text(size=8, angle=30, hjust=1),
axis.title.y=element_blank(), 
plot.title=element_text(size=8, face="bold", hjust=0.5),
axis.title.x=element_blank())+
guides(fill="none")
chloroplot3


fig3 <- plot_grid(chloroplot1,NULL,chloroplot2,chloroplot3, nrow=1, align="h", rel_widths=c(1.1,0.06,1.1,0.9), labels=c("A","","B","C"), label_size=9, vjust=3.7)+
draw_label(0.53, 0.97, label="Symbiont and chlorophyll change", fontface="bold", size=8)+
draw_label(0.33, 0.88, label="Community values\n(inner)", size=6, hjust=1)+
draw_line(x=c(0.35,0.35), y=c(0.72,0.88), size=0.1)+
draw_line(x=c(0.345,0.35), y=c(0.72,0.72), size=0.1)+
draw_line(x=c(0.335,0.35), y=c(0.88,0.88), size=0.1)+
draw_label(label="Amu", 0.28, 0.47, size=5, colour="slategrey", fontface="bold")+
draw_label(label="Asu", 0.3, 0.51, size=5, colour="slategrey", fontface="bold")+
draw_label(label="Gon", 0.3, 0.55, size=5, colour="slategrey", fontface="bold")+
draw_label(label="Ibr", 0.29, 0.585, size=5, colour="slategrey", fontface="bold")+
draw_label(label="Pda", 0.29, 0.63, size=5, colour="slategrey", fontface="bold")+
draw_label(label="Pru", 0.3, 0.765, size=5, colour="slategrey", fontface="bold")+
draw_label(label="Pcy", 0.3, 0.795, size=5, colour="slategrey", fontface="bold")
fig3



