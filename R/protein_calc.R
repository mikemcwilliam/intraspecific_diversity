
filenames <- list.files("data/raw_data/protein_data/txt_files", pattern="*.txt", full.names=TRUE)
ldf <- lapply(filenames, function(x) read.table(x, fill=T, sep="\t"))
names(ldf) <- gsub(".txt", "", filenames)
names(ldf) <- gsub("data/raw_data/protein_data/txt_files/", "", names(ldf) )


store <- NULL

for(i in 1:length(ldf)) {
	
	pda <- ldf[[i]]
	
	# plate
	plate <- pda[c(5:12),c(3:14)]
	colnames(plate) <- seq(1:12)
	rownames(plate) <- LETTERS[seq( from = 1, to = 8 )]
	plate
	
	# flipped plates 
	if(grepl("flip", names(ldf)[i])){
	plate <- plate[seq(dim(plate)[1],1),seq(dim(plate)[2],1)]
	BSA <- c(2, 1, 0.5, 0.25, 0.125) }
	else{
	plate <- plate
	BSA <- c(1.44, 0.72, 0.36, 0.18, 0.09) }
	
	# standards
	standards <- plate[c(1:5),c(1:3)]
	plate[c(1:5),c(1:3)] <- NA
	
	# model
	mod <- lm(rowMeans(standards)~BSA)
	r <- round(summary(mod)$r.squared, 3)
	m <- round(coef(mod)[2], 4)
	b <- round(coef(mod)[1], 3)
	
	# plot
	png(paste("data/raw_data/protein_data/BSA_calibration/", names(ldf)[i], ".png", sep=""))
	plot(x=BSA, y=rowMeans(standards), xlab="Proteins (mg/ml)",
	ylab="Absorbance (%)", main=paste("tray", names(ldf)[i]))
	text(x=0.4, y=0.9, label=substitute(atop(y== m*x + b,R^2 == r), list(m=m, b=b,r=r)))
	lines(x=BSA, y=predict(mod))
	dev.off()
	
	# samples
	samples <- plate
	names(samples) <- rep(c(1:3), 4)
	samples <- rbind(samples[,c(1:3)], 
	samples[,c(4:6)],samples[,c(7:9)], samples[,c(10:12)])
	samples <- samples[c(6:nrow(samples)),]
	samples <- cbind(samples, n=c(1:27), meanOD=rowMeans(samples), m, b, r)
	samples$Tray <- names(ldf)[i]
	samples$Protein <- (samples$meanOD-b)/m   # x=(y-b)/m
	samples$well <- paste(samples$Tray, samples$n, sep=".")
		
	store <- rbind(store, samples)
	}

head(store)

# find mean of a and b
store_b <- store[grepl("b",store$Tray),]
store_b$Tray <- gsub("b", "", store_b$Tray)
store_b$well <- paste(store_b$Tray, store_b$n, sep=".")

store <- store[!grepl("b",store$Tray),]
store$Protein_b <- store_b$Protein[match(store$well, store_b$well)]
store$rsq_b <- store_b$r[match(store$well, store_b$well)]
store$Protein_mgml <- rowMeans(store[,c("Protein", "Protein_b")], na.rm=T)

wells <- read.csv("data/raw_data/protein_data/protein_wells.csv")
wells$well <- paste(wells$Tray, wells$Sample, sep=".")

store$well <- gsub("_flip", "", store$well )

data <- cbind(wells, store[match(wells$well, store$well),])

#write.csv(data, "data/raw_data/protein_data/protein_data.csv")


