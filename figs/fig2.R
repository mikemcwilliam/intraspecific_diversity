
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


esize1 <- ggplot(data=params, aes(x=e.size, y=reorder(name, e.size)))+
geom_vline(xintercept=0, col="black", size=0.1)+
geom_segment(aes(x=e.size+e.se, xend=e.size-e.se, yend=reorder(name, e.size)), size=0.3)+
geom_point(col="white",size=4)+
geom_point(aes(col=type), alpha=0.5, stroke=0.45, size=4)+
geom_point(aes(col=type), shape=21, stroke=0.45, size=4)+
geom_text(aes(label=label), size=1.8, fontface="bold")+
theme_classic()+
ggtitle("Trait ~ Location")+
scale_colour_manual(values=pal2[3:1])+
scale_fill_manual(values=pal2[3:1])+
labs(x="Effect size", y="Trait")+
geom_text(data=NULL, aes(y="Chlorophyll content", x=-2, label=paste("n = ", nrow(tdat), sep="")), size=2.5)+
#xlim(c(-2,8))+
theme(axis.line=element_line(size=0.1), 
#axis.line.y=element_blank(),
#axis.line.x=element_blank(),
plot.title=element_text(size=8, hjust=0.5, face="bold"),
#axis.ticks.x=element_blank(),
plot.background=element_blank(),
axis.text.y=element_text(size=8),
axis.text.x=element_text(size=9, angle=45, hjust=1),
#axis.ticks.y=element_blank(),
axis.title.y=element_text(size=8),
#axis.title.x=element_text(size=8, margin=margin(10,0,0,0))
axis.title.x=element_blank())+
coord_flip()
esize1

nrow(tdat)

s.sizes <- data.frame(table(tdat$Species))
params2$ss <- s.sizes$Freq[match(params2$sp, s.sizes$Var1)]

sum(s.sizes$Freq)

params3 <- params2
params3$sp[params3$sp=="Goniastrea spp."] <- "Goniastrea sp. "
params3$e.seX <- ifelse(params3$esize <0, params3$e.se*-1, params3$e.se)
params3$ss1 <- paste("n = ", params3$ss, sep="")

esize.sp<-ggplot(params3, aes(label, esize))+
geom_bar(stat="identity", aes(col=type,fill=type), size=0.2, width=0.8)+
geom_bar(data=params3[params3$sig=="*",],aes(fill=type), stat="identity",size=0.25, width=0.8, col="black")+
geom_segment(aes(y=esize, yend=esize+e.seX, x=label, xend=label, col=type), size=0.3)+
#geom_text(aes(label=sig), vjust=0.75)+
facet_wrap(~sp, nrow=2, scale="free_x")+
geom_text(data=unique(params3[,c("ss1", "sp")]), aes(x="ZD", y=-4, label=ss1), size=2.5)+
coord_flip()+
guides(col="none")+
#ggtitle("Species-specific responses")+
#scale_y_continuous(breaks=c(-6,-4,-2,0,2,4,6,8))+
geom_hline(yintercept=0, size=0.1)+
scale_fill_manual(values=pal2[3:1])+
scale_colour_manual(values=pal2[3:1])+
scale_y_continuous(breaks=c(-4, -2,0,2, 4),limits=c(-5.5,5.5))+
theme_classic()+guides(fill="none")+theme(strip.background=element_blank(), axis.text.x=element_text(size=8), 
axis.text.y=element_text(size=9), 
axis.line=element_line(size=0.2), 
axis.title=element_text(size=8),
strip.text=element_text(size=8, face="italic"),
plot.title=element_text(size=8, hjust=0.5, face="bold", margin=margin(3,0,12,0)),
axis.ticks=element_line(size=0.2))+
labs(y="Effect size", x="Trait")
esize.sp

sp.plot<-plot_grid(NULL, esize.sp, ncol=1, rel_heights=c(0.05,1))+
draw_grob(stag, 0.06, 0.6, width=0.13, height=0.08)+
draw_grob(card, 0.28, 0.6, width=0.13, height=0.08)+
draw_grob(gon, 0.51, 0.6, width=0.13, height=0.09)+
draw_grob(brug, 0.74, 0.58, width=0.13, height=0.1)+
draw_grob(pdam, 0.06, 0.15, width=0.13, height=0.08)+
draw_grob(cyl, 0.29, 0.15, width=0.13, height=0.08)+
draw_grob(rus, 0.52, 0.13, width=0.13, height=0.1)+
draw_label(label="Species trait responses", 0.5, 0.95, fontface="bold", size=8)
sp.plot

legend2 <-get_legend(esize1+theme(legend.background=element_blank(),
legend.text=element_text(size=8),
#legend.position=c(0.2,0.9),
legend.position=c(0,0.3),
legend.title=element_blank(),
legend.key.size=unit(0.7,"mm"))+
guides(col = guide_legend(direction="vertical", override.aes = list(size = 2)))
)

fig.2 <- plot_grid(esize1+guides(col="none"), sp.plot, NULL, legend2, NULL, nrow=1,
rel_widths=c(0.38,1, -0.1, 0.1,0.05), labels=c("A","B"),label_size=9)
fig.2

