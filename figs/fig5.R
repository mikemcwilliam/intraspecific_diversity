


ord <- c("Zooxanthellae density","Chlorophyll content","Tissue biomass","Protein biomass","Skeletal density","Surface:volume ratio","Branch density","Branch width","Corallite size")

SSlong$variableX <- ifelse(SSlong$variable=="esize","Effect size",ifelse(SSlong$variable=="Prop_intra","Intraspecific", ifelse(SSlong$variable =="Prop_turn","Species composition",ifelse(SSlong$variable=="Prop_cov","Covariation", NA))))
SSlong$variableX <- factor(SSlong$variableX, levels=c("Effect size","Intraspecific","Species composition","Covariation"))

datfig5 <- subset(SSlong[SSlong$SS=="SStot",], variableX!="Effect size")
datfig5$name <- info$name[match(datfig5$label, info$label)]
datfig5$name <- factor(datfig5$name, levels=ord)

decomp <- ggplot(data=datfig5)+geom_bar(stat="Identity", aes(x=value, y=name, fill=type))+
geom_vline(xintercept=0, size=0.1)+
facet_wrap(~variableX, scales="free_x")+
scale_fill_manual(values=pal2[3:1])+
labs(x="Variance explained")+
#guides(fill="none")+
#scale_x_continuous(expand=c(0,0))+
theme_classic()+
theme(axis.line.x=element_line(size=0.1), 
axis.line.y=element_blank(),
axis.title.y=element_blank(), 
legend.position="top",
#legend.position=c(0.7,0.3),
legend.title=element_blank(),
legend.background=element_blank(),
legend.margin=margin(1,0,-7,0),
strip.text=element_text(size=8),
strip.background=element_blank(),
legend.key.height=unit(0.5,"mm"),
plot.background=element_blank(),
legend.key.width=unit(1,"mm"),
legend.text=element_text(size=8),
plot.title=element_text(hjust=0.5, size=8, face="bold"),
axis.text=element_text(size=8), 
axis.title.x=element_text(size=8))
decomp 



pdat2 <- melt(pdat)

pdat2$variable2 <- ifelse(pdat2$variable=="fixed", "Species composition", "Intraspecific")
pdat2$variable2<-factor(pdat2$variable2, levels=c("Species composition", "Intraspecific"))
pdat2$type2 <- ifelse(pdat2$variable=="intra", "intra", as.character(pdat2$type))


legendplot <- ggplot()+
geom_bar(data=pdat2[pdat2$variable %in% c("intra","fixed"),], aes(value, reorder(name, -value), fill=variable2), stat="identity", size=0.1,  col="black")+
theme_classic()+
#scale_fill_manual(values=c(pal2[3], "white"))+
scale_fill_manual(values=c("grey", "white"))+
theme(legend.position="right", 
legend.title=element_blank(), 
legend.key.width=unit(3,"mm"),
legend.key.height=unit(3,"mm"),
legend.text=element_text(size=8))

dat5 <- pdat2[pdat2$variable %in% c("intra","fixed"),]
dat5$name <- factor(dat5$name, levels=ord)

cwm.es <- ggplot(data=dat5)+
#annotation_custom(get_legend(legendplot), xmin=3.9,xmax=4.4, ymin=8, ymax=10)+
geom_bar(data=dat5, aes(value, name, fill=type2, col=type), stat="identity", size=0.5,  width=0.8)+
geom_vline(xintercept=0, size=0.1)+
geom_point(data=pdat, aes(specific, reorder(name, -specific)), shape=8, size=0.75, stroke=0.25)+
scale_fill_manual(values=c("white",pal2[3:1]))+
scale_colour_manual(values=pal2[3:1])+
guides(fill="none", col="none")+
labs(x="Effect size")+
theme_classic()+
theme(axis.line=element_line(size=0.1), 
axis.title.y=element_blank(), 
plot.title=element_text(hjust=0.5, size=8),
axis.text.y=element_blank(),
plot.background=element_blank(),
axis.text.x=element_text(size=8), 
axis.title.x=element_text(size=8))
cwm.es

fig.5 <- plot_grid(NULL,
plot_grid(NULL, get_legend(decomp),rel_widths=c(-0.1,1)),
NULL,
plot_grid(
plot_grid(decomp+guides(fill="none")+ggtitle(""), cwm.es+ggtitle("CWM ~ location"), nrow=1, rel_widths=c(1, 0.35, 0.1), labels=c("A","B"), label_size=9, align="h", axis="bt"),NULL,
get_legend(legendplot), rel_widths=c(1,-0.11,0.2), nrow=1),
ncol=1, rel_heights=c(0.1,0.1,-0.1,1))+
draw_label(label="Source of trait variation", x=0.5, y=0.95, fontface="bold", size=8)
fig.5


