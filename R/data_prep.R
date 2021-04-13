

# Coral Traits Database means  - www.coraltraits.org
#source("data/traitdatabase/ctb_prep.R")
ctb<-read.csv("data/traitdatabase/output.csv", as.is=TRUE)
ctb<-ctb[ctb$Zooxanthellate=="zooxanthellate",] # zooxanthellate only
ctb<-ctb[!is.na(ctb$Abundance.GBR),] # gbr only
nrow(ctb)

# trait biogeography - McWilliam et al. 2018 PNAS
dat<-read.csv("data/pnas2018/traitbiogeography.csv", as.is=TRUE)
dat$species<-gsub("Coscinarea", "Coscinaraea",dat$species)
dat$species<-gsub("Dipsastrea", "Dipsastraea",dat$species)
dat$species<-gsub("erythrea","erythraea",dat$species)
dat$species<-gsub("Stylarea","Stylaraea",dat$species)
dat$species<-gsub("Paramontastrea", "Paramontastraea",dat$species)
dat$species<-gsub("Montastrea", "Montastraea",dat$species)
#dat$species<-gsub("Helioseris cuc", "Leptoseris cuc",dat$species)
dat$species<-gsub("Caulastrea", "Caulastraea",dat$species)
dat$species<-gsub("Anomastrea", "Anomastraea",dat$species)
dat$species<-gsub("Favia leptophylla", "Mussismilia leptophylla",dat$species)

# TRAIT MISTAKES..
# Tables should be given highest SA & space (original 2, 2, 2)
dat$cat_SA_vol[dat$raw_growth_form=="tables_or_plates"]<-5
dat$cat_spacesize[dat$raw_growth_form=="tables_or_plates"]<-4
dat$cat_colonyheight[dat$raw_growth_form=="tables_or_plates"]<-3
# branching closed? Bushy or staghorn?  (original 4, 3, 4)
#dat$cat_spacesize[dat$raw_growth_form=="branching_closed"]<-2
#dat$cat_SA_vol[dat$raw_growth_form=="branching_closed"]<-4
#dat$cat_colonyheight[dat$raw_growth_form=="branching_closed"]<-3

# match pnas and ctb
matches<-dat$species[match(ctb$species, dat$species)]
ctb$species[is.na(matches)] # ctb species NOT found in pnas data
# find out where these 4 are.. but remove for now
ctb<-ctb[!is.na(matches),]
nrow(ctb)
dat<-dat[match(ctb$species, dat$species),]
nrow(dat)

# select traits of interest
colnames(ctb)
raw<-ctb[,c("species", "Genus.fossil.age","Species.age.phylogeny","Abundance.GBR", "Range.size", "Growth.form.typical","Corallite.width.maximum","Colony.maximum.diameter", "Growth.rate","Oocyte.size.at.maturity","Skeletal.density","Substrate.attachment", "Coloniality", "Mode.of.larval.development","Sexual.system", "Symbiodinium.sp..in.propagules", "Water.clarity.preference", "Polyps.per.area")]
head(raw)

colnames(dat)
cats<-dat[,c("species", "genus","domain","cat_growthrate", "cat_corallitesize", "cat_colonydiameter" , "cat_skeletaldensity", "cat_colonyheight","cat_SA_vol", "cat_spacesize", "dat_growth", "dat_skeletal", "dat_corallite" ,"dat_colonydiameter", "reproductive_mode")]
head(cats)
head(raw)

dat<-merge(raw, cats, by="species")
head(dat)
tail(dat)
nrow(dat)

# clipperton species
clipperton<-c("Porites lobata", "Pavona varians", "Pavona minuta", "Leptoseris scabra", "Pavona maldivensis", "Pocillopora meandrina", "Porites lutea")
dat$clipperton<-NA
dat$clipperton<-ifelse(dat$species %in% clipperton, 1, 0)

# bleaching susceptibility... 
bri<-read.csv("data/BRI/bleaching_response_index.csv")[,c(1,2)]
names(bri)<-c("species", "BRI")
dat$BRI_genus <- bri$BRI[match(dat$genus, bri$species)]
dat$BRI <- bri$BRI[match(dat$species, bri$species)]
dat$BRI[is.na(dat$BRI)] <- dat$BRI_genus[is.na(dat$BRI)]

#head(dat)
#hist(log10(dat$BRI))

# LAtitude

lat <- read.csv("data/data_20210211.csv", as.is=TRUE)[c("specie_name", "trait_name", "value")]
lat1 <- lat[lat$trait_name=="Northern-most range edge", c("specie_name", "value")]
names(lat1) <- c("species", "north")
lat2 <- lat[lat$trait_name=="Southern-most range edge", c("specie_name", "value")]
names(lat2) <- c("species", "south")

lat <- merge(lat1, lat2)
lat$lat_range <- lat$north - lat$south

dat <- merge(dat, lat, all.x=TRUE)

# restorability from PLoS paper
dat$restore <- 0
dat$restore[dat$Growth.form.typical %in% c("branching_open", "branching_closed")] <- 6
dat$restore[dat$Growth.form.typical %in% c("massive", "submassive")] <- 5
dat$restore[dat$Growth.form.typical %in% c("laminar")] <- 4
dat$restore[dat$Growth.form.typical %in% c("encrusting", "encrusting_long_uprights")] <- 3
dat$restore[dat$Growth.form.typical %in% c("digitate", "corymbose", "tables_or_plates")] <- 2
dat$restore[dat$Growth.form.typical %in% c("columnar", "hispidose")] <- 1
dat$restore <- (dat$restore) / 6

# dat$restore[dat$Growth.form.typical %in% c("branching_open", "massive", "encrusting", "branching_closed", "tables_or_plates", "encrusting_long_uprights")] <- 1
# dat$restore[dat$Growth.form.typical %in% c("digitate", "corymbose", "laminar")] <- 0.5
# dat$restore[dat$Growth.form.typical %in% c("columnar", "hispidose", "submassive")] <- 0.1

# unique(dat$genus)
# dat$restore2 <- 1
# dat$restore2[dat$genus %in% c("Acropora")] <- 30


# Families
fam <- read.csv("data/species-5.csv", as.is=TRUE)[c("master_species", "family_molecules", "family_morphology")]
names(fam) <- c("species", "family_molecules", "family_morphology")
dat <- merge(dat, fam, all.x=TRUE)
unique(dat$family_morphology)


# Normalise
dat$range <- dat$Range.size / max(dat$Range.size)
dat$abund[dat$Abundance.GBR=="common"] <- 1
dat$abund[dat$Abundance.GBR=="uncommon"] <- 0.5
dat$abund[dat$Abundance.GBR=="rare"] <- 0.1
dat$bleach <- 1 - (dat$BRI / 100)

dat$lat_range <- dat$lat_range / max(dat$lat_range)
# dat$lat_poleward <- apply(dat[c("north", "south")], 1, function(x) max(abs(x)))
dat$lat_poleward <- abs(dat$south) / max(abs(dat$south))

# Anonymous species

dat <- dat[c("species", "range", "lat_poleward", "abund", "bleach", "restore", "cat_growthrate", "cat_corallitesize", "cat_colonydiameter" , "cat_skeletaldensity", "cat_colonyheight","cat_SA_vol", "cat_spacesize", "clipperton")]
dat <- na.omit(dat)
dim(dat)
dat$species_n <- paste0("Species ", 1:nrow(dat))

dat[c("species", "clipperton")]


# export
write.csv(dat, "data/data.csv")

