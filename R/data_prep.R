

####################################
####################################
####################################
# ENVIRONMENTS

env <- subset(read.csv("data/raw_data/sites.csv"), select=-c(Notes))
env$Date <- as.Date(gsub("_", "-", env$Date))

# depth 1m only
env <- subset(env, Depth_m==1)

# create site ID
env$SiteID <- ifelse(env$Site=="Venture Lake", "H", ifelse(env$Site=="Venture Lake outside", "G",
ifelse(env$Site=="Nikko Bay", "F",
ifelse(env$Site=="Mutmelachel", "D",
ifelse(env$Site=="Falcon Reef", "E",
ifelse(env$Site=="Pats Dropoff", "C",
ifelse(env$Site=="Eastern Barrier", "B",
ifelse(env$Site=="Ulong Rock", "A", NA))))))))

# new longitude for mapping
env$long2 <-ifelse(env$SiteID=="H",134.481, env$long)

# write.csv(env, "data/envs.csv")

####################################
####################################
####################################
# ABUNDANCE

lit <- read.csv("data/raw_data/abundance/PalauSurveys.csv")
head(lit)

# sum colonies
lit$cover <- rowSums(lit[,grepl("X", colnames(lit))],na.rm=T)/1000*100

# n colonies
n <- ifelse(lit[,grepl("X", colnames(lit))] > 0, 1, 0)
lit$n <- rowSums(n, na.rm=T)

# remove zeros
sp.total <- aggregate(cover~Species, lit, sum)
zeros <- sp.total$Species[sp.total$cover==0]
lit <- lit[!lit$Species %in% zeros,]

# remove intercept data
lit <- lit[,!grepl("X", colnames(lit))]

# add site ID
lit$SiteID <- env$SiteID[match(lit$Site, env$Site)]

# write.csv(lit, "data/abundance.csv")

####################################
####################################
####################################
# TRAITS 

dat <- read.csv("data/raw_data/samples.csv")
head(dat)
#dat <- dat[!dat$Photo.2=="IMG_1759",] # deep

####################################
### Fragment area/volume/density

sav <- read.csv("data/raw_data/3d_files/3d_data.csv")

# average duplicates
sav$ID[duplicated(sav$ID)]
sav <- aggregate(cbind(full_VOLmm3,full_SAmm2, live_SAmm2, Mass_g) ~ ID, sav, mean)

# align/attach
dat[,names(sav)] <- sav[match(dat$ID, sav$ID), names(sav)]
dat[is.na(dat$live_SAmm2),"ID"] # missing data

# Skeletal Density
dat$SD_gcm3 <- dat$Mass_g / (dat$full_VOLmm3 / 1000)

# Surface area to volume
dat$SAV <- dat$full_SAmm2/dat$full_VOLmm3

####################################
### Proteins 

# source("R/protein_calc.R")

pro <- read.csv("data/raw_data/protein_data/protein_data.csv") 
head(pro)

# convert new 660 readings
pro$P_conc <- ifelse(pro$Dye == "NEW_660", pro$Protein_mgml*10, pro$Protein_mgml)

# average duplicates (within dyes) - only s7_prus_1
pro <- aggregate(P_conc ~ ID + Dye, pro, mean)

# align/attach - new (+ old)
dat$Protein_mgml <- pro$P_conc[match(dat$ID, pro$ID)] 
dat[is.na(dat$Protein_mgml) ,"ID"] # missing data

####################################
### Chlorophyll 

# source("data/chlorophyll_data/chlorophyll_calc.R")

chl <- read.csv("data/raw_data/chlorophyll_data/chlorophyll_data.csv") 
head(chl)

# duplicates
chl$ID[duplicated(chl$ID)]
chl <-aggregate(.~ID, chl[,c("ID","chl_a","chl_c")], mean)

# align/attach
dat[,c("chl_a","chl_c")] <- chl[match(dat$ID, chl$ID),c("chl_a","chl_c")]
dat[is.na(dat$chl_a),"ID"] # missing data

####################################
### Zooxanthellae 

# attach IDs to images
ref<-read.csv("data/raw_data/zoox_data/image_reference.csv")

zoo<-read.csv("data/raw_data/zoox_data/zoox_counts.csv")
zoo$ID<-ref$ID[match(zoo$Slice, ref$Image)]
head(zoo)

zoo$zoox <- as.numeric(as.character(zoo$zoox))

# labels/duplicates
data.frame(table(zoo$ID))
subset(data.frame(table(zoo$ID)), Freq > 6)
zoox <- aggregate(cbind(Auto_Count, Average.Size, zoox) ~ ID, zoo, mean)

# align/attach
#dat$zoox_old<-zoox$Auto_Count[match(dat$ID, zoox$ID)]
dat$zoox<-zoox$zoox[match(dat$ID, zoox$ID)]
dat[is.na(dat$zoox),"ID"] # missing data 

####################################
### Biomass

bio<-read.csv("data/raw_data/biomass/biomass.csv")
head(bio)
#bio[grep("excess",bio$ID),]

# duplicates
bio$ID[duplicated(bio$ID)]
bio <- aggregate(cbind(AFDW_5ml, Ash.Weight) ~ ID + top, bio, mean)

# align/attach
dat[,"AFDW_5ml"] <- bio[match(dat$ID, bio$ID),"AFDW_5ml"] 
dat[is.na(dat$AFDW_5ml) ,"ID"] # missing data

####################################
# Colours 

dist <- read.csv("data/raw_data/colours/dist.csv", row.names=1)
#pca <- cmdscale(as.dist(dist))

col <- read.csv("data/raw_data/colours/coldat.csv")
#df[,c(1:3)] <- round(df[,c(1:3)]*255, 0)

isIDmax <- with(col, ave(rgb.Pct, IMG, FUN=function(x) seq_along(x)==which.max(x)))==1
maxcol <- col[isIDmax, ]

dat$col <- maxcol$col[match(dat$Photo.1, maxcol$IMG)]
dat[is.na(dat$col) ,"ID"] # missing data

dat$col <- dat$col
dat$col[dat$ID=="s1_anac_2"] <- dat$col[dat$ID=="s1_anac_1"]
dat$col[dat$ID=="s2_card_1"] <- dat$col[dat$ID=="s3_card_2"]
dat$col[dat$ID=="s3_card_4"] <- dat$col[dat$ID=="s3_card_3"]

####################################
### Morphology from images

mor <- read.csv("data/raw_data/morph_photos/morphology.csv")

scl <- 2.5   # 1 scalebar box = 2.5 cm

mor$colony_area <- mor$planar_area * (scl^2)
mor$branch_density <- mor$n_branches / (mor$n_in_area * (scl^2))
mor$branch_width <- rowMeans(mor[,grepl("branch_width", colnames(mor))], na.rm=T) * scl
mor$branch_height <- rowMeans(mor[,grepl("branch_height", colnames(mor))], na.rm=T) * scl
mor$polyp_size <- mor$corallite_width * scl
mor$polyp_density <- mor$n_corallites / ((mor$box_width^2)*(scl^2))
mor$base_ratio <- mor$top_width / mor$base_width 
mor$height <- mor$colony_height * scl

morphs <- c("colony_area","branch_density", "branch_width", "branch_height", "polyp_size", "polyp_density", "base_ratio", "height")
dat[,morphs] <- mor[match(dat$Photo.1, mor$IMG),morphs]

####################################
####################################
####################################
### NORMALISE/SCALE

dat$chl_a <- dat$chl_a*(9/5)*15 # ul/samp
dat$chl_c <- dat$chl_c*(9/5)*15 # ul/samp
dat$chl_t <- dat$chl_a + dat$chl_c
dat$zoox <- (dat$zoox*(10^4)*15)/10^6 # cells/samp
dat$protein <- dat$Protein_mgml*15 #mg/samp
dat$biomass <- dat$AFDW_5ml*1000 *3 #mg/samp

#names(dat)
dat_norm <- dat
norms <- c("protein","chl_a", "chl_c","chl_t", "zoox", "biomass")
dat_norm[,norms] <- dat_norm[,norms] / ((dat$live_SAmm2)/100)

# write.csv(dat_norm, "data/traits.csv")

####################################
####################################
####################################
### PORITES MASSIVE

por <- dat_norm[1,] 
por[,2:ncol(por)] <- NA
por[,c("Site.N", "Site", "Location","Bag","Photo.1", "Photo.2", "Tag", "ID")] <- "CTB"
por[,"Genus"] <- "Porites"
por[,"Species"] <- "Porites lutea"

pvals <- read.csv("data/raw_data/massive_porites/massive_porites.csv")
pmean <- aggregate(Value~Trait, pvals, mean)
por$zoox <- pmean$Value[pmean$Trait=="zoox_den"]
por$protein <- pmean$Value[pmean$Trait=="protein_mass"]
por$protein_av <- pmean$Value[pmean$Trait=="protein_mass"]
por$biomass <- pmean$Value[pmean$Trait=="total_biomass"]
por$chl_a <- pmean$Value[pmean$Trait=="chl_a"]
por$chl_c <- pmean$Value[pmean$Trait=="chl_a"]
por$chl_t <- pmean$Value[pmean$Trait=="chl_a"]
por$SD_gcm3 <- pmean$Value[pmean$Trait=="skel_den"]

lobfill <- c("branch_density", "branch_width","branch_height", 'SAV')
lobs <- dat_norm[dat_norm$Genus=="Lobophyllia",]
lobs <- aggregate(.~Genus, lobs[,c("Genus", lobfill)], mean)
por[,lobfill] <- lobs[,lobfill]
por

# write.csv(por, "data/massive_porites/por.csv")


####################################
####################################
####################################
# SPECIES

spp <- read.csv("data/raw_data/species.csv")
n.spp<-data.frame(table(dat$Species))
spp$samples<-n.spp$Freq[match(spp$species, n.spp$Var1)]
occdat<-subset(data.frame(table(dat$Species, dat$Site)), Freq>0)
occdat<-data.frame(table(occdat$Var1))
spp$occurence<-occdat$Freq[match(spp$species, occdat$Var1)]
dat$occurence<-occdat$Freq[match(dat$Species, occdat$Var1)]
spp$Taxa <- dat$Taxa[match(spp$species, dat$Species)]
spp$samp_site <- spp$samples/spp$occurence

locs<-as.data.frame.matrix(table(dat$Species, dat$Location))
locs$Occ <- ifelse(locs[,1]==0, "Outer", ifelse(locs[,3]==0, "Inner", "Both"))
spp[,colnames(locs)] <- locs[match(spp$species, rownames(locs)),]
spp

# write.csv(spp, "data/species.csv")


