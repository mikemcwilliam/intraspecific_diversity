
library("rgdal")
library("ggsn") #ggplot scalebars
library("rworldmap")

# Palau shapefiles

shp<-readOGR("data/raw_data/shp_files/pw_unepwcmc_all_coralreefs.shp")
unique(shp$rb_attrib)

land <- fortify(shp[shp$rb_attrib == "land",])
br <- fortify(shp[shp$rb_attrib == "barrier island",])
fr <- fortify(shp[shp$rb_attrib == "fringing island",])
pr <- fortify(shp[shp$rb_attrib == "patch island",])
reef <- rbind(br, fr, pr)

########################################
# worldmap

worldMap <- getMap()
world.points <- fortify(worldMap)
world.points$region <- world.points$id
world.df <- world.points[,c("long","lat","group", "region")]

x_lines <- seq(-180,180, by = 20)
y_lines <- seq(-80,80, by = 20)

worldmap <- ggplot() + 
	   geom_segment(aes(y = -160, yend = 180, x = x_lines, xend = x_lines), linetype = "solid", col="black", size=0.05) +
	    geom_segment(aes(x = -100, xend = 205, y = y_lines, yend = y_lines), linetype = "solid", col="black", size=0.05) +
  geom_polygon(data = world.df, aes(x = long, y = lat, group = group), fill="white", col="black", size=0.05) +
  geom_point(data=NULL, aes(x=134.4772, y=7.3411), col="black", shape=22)+
  coord_map("ortho", orientation=c(20, 125, 0))+
  theme_void()
worldmap  

########################################
# NOAA SST data
  
# https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMBsstdmday
#turl <- "https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMBsstdmday.csvp?sst%5B(2006-01-16T12:00:00Z):1:(2020-12-16T12:00:00Z)%5D%5B(0.0):1:(0.0)%5D%5B(6.8):1:(8.2)%5D%5B(134):1:(135)%5D"  
#download.file(turl, "data/raw_data/SST.csv", mode='wb')

sst <- read.csv("data/raw_data/SST.csv")
names(sst) <- c("time","alt", "lat", "long", "sst")

head(sst)

sst$time <- as.Date(sst$time)
sst$month <- substr(sst$time, 6,7)
sst$year <- substr(sst$time, 1,4)
head(sst)

seasons <-data.frame(x=as.Date(c("2011-02-01","2011-12-01","2012-12-01","2013-12-01","2014-12-01","2015-12-01","2016-12-01","2017-12-01","2018-12-01","2019-12-01","2020-12-01")), 
xend=as.Date(c("2011-04-01","2012-04-01","2013-04-01","2014-04-01","2015-04-01","2016-04-01","2017-04-01","2018-04-01","2019-04-01","2020-04-01","2021-01-01")),
y=Inf)

mid <- c(7.3, 134.383) # midpoint
a <- 0.18 # area 
 
sst <- subset(sst, lat>mid[1]-a/2-0.05 & lat<mid[1]+a/2+0.05)
sst <- subset(sst, long>mid[2]-a-0.2 & long<mid[2]+a+0.2)

sstloc <- rbind(
cbind(subset(sst, lat>7.3 & lat<7.35 & long>134.45 & long <134.51), loc="Inner"),
cbind(subset(sst, lat>7.349 & lat<7.351 & long>134.5 & long <134.54), loc="Inner"),
cbind(subset(sst, lat>7.25 & lat<7.3 & long>134.49 & long <134.55), loc="Outer"),
cbind(subset(sst, lat>7.3 & lat<7.35 & long>134.3 & long <134.4), loc="Lagoon"))
head(sstloc)
tempdat<-subset(aggregate(sst~time+loc+year,sstloc, mean), year>2010)

pal <- wes_palette("Zissou1", 3, type = "continuous")

ggplot()+
geom_rect(data=NULL, aes(xmin=min(seasons$x),xmax=max(seasons$x), ymin=-Inf, ymax=Inf), alpha=0.6, fill="white")+
geom_rect(data=seasons, aes(xmin=x, xmax=xend, ymin=y, ymax=-y), fill="white", alpha=1)+
geom_line(data=tempdat, aes(time, sst, col=loc), size=0.18)

##################################
# Map

sub <- sst[sst$year > 2014 & sst$sst < 32,]
df <- aggregate(sst~lat+long, sub, mean, na.rm=T)
head(df)

ewbrks <- seq(134.2,134.5, 0.1)
nsbrks <- seq(7.25,7.35,0.05)
ewlbls <- unlist(lapply(ewbrks, function(x) ifelse(x < 0, paste(x, "°E"), ifelse(x > 0, paste(x, "°E"),x))))
nslbls <- unlist(lapply(nsbrks, function(x) ifelse(x < 0, paste(x, "°S"), ifelse(x > 0, paste(x, "°N"),x))))

map <- ggplot( )+ 
  geom_raster(data = df, aes(x =long, y = lat, fill = sst), interpolate=T, alpha=1) +
  geom_polygon(data=reef, aes(long,lat, group=group), fill="#d0dede" , size=0.05, col="black")+
  geom_polygon(data=land, aes(long,lat, group=group), fill="grey65", size=0.1, col="black")+
  scale_fill_distiller(palette="YlGn", direction=1, breaks=c(29.5, 30), guide=guide_colourbar(ticks.colour = "black", frame.colour="black", title.vjust=-0.5,ticks.linewidth = 1))+
  coord_cartesian(c(mid[2]-a, mid[2]+a) ,c(mid[1]-a/2, mid[1]+a/2))
  map