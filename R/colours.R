
library("colordistance")

rgblist <- getHistList("data/colours/colours_img")
rgbdat<- do.call(rbind, rgblist)
rgbdat$IMG <- substr(rownames(rgbdat), 1,8)
rownames(rgbdat) <- c(1: nrow(rgbdat))
rgbdat$col <- apply(rgbdat, 1, function(x){rgb(x[1], x[2], x[3], maxColorValue=1)})

head(rgbdat)

hsvlist <- getHistList("data/colours/colours_img", hsv=TRUE)
hsvdat<- do.call(rbind, hsvlist)
hsvdat$IMG <- substr(rownames(hsvdat), 1,8)
rownames(hsvdat) <- c(1: nrow(hsvdat))

head(hsvdat)

coldat <- cbind(rgbdat[,c("IMG", "col")], rgbdat[,c("r","g","b")], rgb.Pct = rgbdat$Pct, hsvdat[,c("h","s","v")], hsv.Pct = hsvdat$Pct)
head(coldat)

#write.csv(coldat, "data/data_files/colours/coldat.csv")


# EMD distance between all image pairs

dist <- imageClusterPipeline("data/colours/colours_img", color.space = "rgb", distance.method = "emd", cluster.method = "hist", plot.heatmap = FALSE)


#write.csv(dist, "data/colours/dist.csv")



