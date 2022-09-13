

filenames <- list.files("data/chlorophyll_data", pattern="*.txt", full.names=TRUE)

ldf <- lapply(filenames, function(x) read.table(x, fill=T, sep="\t"))
names(ldf) <- gsub(".txt", "", filenames)
names(ldf) <- gsub("data/chlorophyll_data/plate", "", names(ldf) )

chl<-read.csv("data/chlorophyll_data/chlorophyll_wells.csv")
chl$well <- paste(gsub("Plate ", "", chl$Tray), chl$Sample, sep=".")

store <- NULL

for(i in 1:length(ldf)) {	
	pda <- ldf[[i]]	
	
	# 3 wavelength plates
	plates <- list(pda[c(5:12),c(3:14)], 
	pda[c(5:12),c(16:27)], pda[c(5:12),c(29:40)])
	names(plates) <- c("A630", "A663", "A750")	
	
	df <- data.frame(well = paste(names(ldf)[i], c(1:32), sep="."))
	df$ID <- chl$ID[match(df$well, chl$well)]
	
	for (n in 1:length(plates)){
		
		# arrange by sample
		samples <- plates[[n]]
		names(samples) <- rep(c(1:3), 4)
		samples <- rbind(samples[,c(1:3)], samples[,c(4:6)],
		samples[,c(7:9)], samples[,c(10:12)])
		
		# find aborbance means
		samples$mean <- rowMeans(samples[,c("1","2","3")])
		
		# find blank mean
		samples$blanks <- mean(samples[df$ID=="Acetone","mean"], na.rm=T)
		
		# difference from blank
		samples$Abs <- samples$mean - samples$blanks
		names(samples) <- paste(names(plates)[n], names(samples), sep="_")
		
		df <- cbind(df, samples[,c(4:6)])
		}
	 store <- na.omit(rbind(store, df))
	} 
	
# calculate chlorophyll 

		ml.ext <- 15 
		ml.acet <- 4 

store$chl_a <- ((11.43*(store$A663_Abs-store$A750_Abs))-(0.64*(store$A630_Abs-store$A750_Abs)))

store$chl_c<-((-3.63*(store$A663_Abs-store$A750_Abs))+(27.09*(store$A630_Abs-store$A750_Abs)))

# find mean of a and b

store_mean <- store
store_mean$well <- gsub("b", "", store_mean$well)
store_mean <- aggregate(.~well+ID, store_mean, mean)
head(store_mean)
			
		store_mean <- store_mean[!store_mean$ID=="Acetone",]
		
		write.csv(store_mean, "data/chlorophyll_data/chlorophyll_data.csv")

