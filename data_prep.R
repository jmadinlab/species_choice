


source("data/pnas2018/pnas_prep.R")
source("data/traitdatabase/ctb_prep.R")


# Coral Traits Database means  - www.coraltraits.org
ctb<-read.csv("data/traitdatabase/output.csv")
ctb<-ctb[ctb$Zooxanthellate=="zooxanthellate",] # zooxanthellate only
ctb<-ctb[!is.na(ctb$Abundance.GBR),] # gbr only
nrow(ctb)


# trait biogeography - McWilliam et al. 2018 PNAS
pnas<-read.csv("data/pnas2018/output.csv")
nrow(pnas)
matches<-pnas$species[match(ctb$species, pnas$species)]
ctb$species[is.na(matches)] # ctb species NOT found in pnas data
# find out where these 4 are.. but remove for now

# narrow ctb data to that that match 
ctb<-ctb[!is.na(matches),]
nrow(ctb)

# narrow down pnas data to gbr species
pnas<-pnas[match(ctb$species, pnas$species),]
nrow(pnas)

# select traits of interest
colnames(ctb)
raw<-ctb[,c("species", "Genus.fossil.age","Species.age.phylogeny","Abundance.GBR", "Range.size", "Growth.form.typical","Corallite.width.maximum","Colony.maximum.diameter", "Growth.rate","Oocyte.size.at.maturity","Skeletal.density","Substrate.attachment", "Coloniality", "Mode.of.larval.development","Sexual.system", "Symbiodinium.sp..in.propagules", "Water.clarity.preference", "Polyps.per.area")]

colnames(pnas)
cats<-pnas[,c("species", "genus","domain","cat_growthrate", "cat_corallitesize", "cat_colonydiameter" , "cat_skeletaldensity", "cat_colonyheight","cat_SA_vol", "cat_spacesize", "dat_growth", "dat_skeletal", "dat_corallite" ,"dat_colonydiameter", "reproductive_mode")]
head(cats)
head(raw)

dat<-merge(raw, cats, by="species")
head(dat)
nrow(dat)


# clipperton species
clipperton<-c("Porites lobata", "Pavona varians", "Pavona minuta", "Leptoseris scabra", "Pavona maldivensis", "Pocillopora meandrina", "Porites lutea")
dat$clipperton<-NA
dat$clipperton<-ifelse(dat$species %in% clipperton, 1, 0)








# export
write.csv(dat, "data/data.csv")

