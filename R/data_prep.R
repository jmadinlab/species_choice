

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
dat$BRI<-bri$BRI[match(dat$species, bri$species)]
#head(dat)
#hist(log10(dat$BRI))

# restorability 
dat$restore <- 0
dat$restore[dat$Growth.form.typical %in% c("branching_open", "massive", "encrusting", "branching_closed", "tables_or_plates", "encrusting_long_uprights")] <- 1
dat$restore[dat$Growth.form.typical %in% c("digitate", "corymbose", "laminar")] <- 0.5
dat$restore[dat$Growth.form.typical %in% c("columnar", "hispidose", "submassive")] <- 0.1

# Families
fam <- read.csv("data/species-5.csv", as.is=TRUE)[c("master_species", "family_molecules", "family_morphology")]
names(fam) <- c("species", "family_molecules", "family_morphology")
dat <- merge(dat, fam, all.x=TRUE)
unique(dat$family_morphology)

# dat$BRI <- dat$BRI/100
# mod <- glm(BRI ~ Growth.form.typical + family_molecules, dat, family=binomial)
# drop1(mod, test="Chisq")
# hist(dat$BRI)

# tax<-read.csv("data/taxonomy.csv")
#head(tax)

# miz <- read.csv("data/mizerek/338_2018_1702_MOESM1_ESM.csv", as.is=TRUE)
# suc <- miz$BI * miz$Number.of.colonies
# fai <- miz$Number.of.colonies - suc
# suc <- round(suc)
# fai <- round(fai)
# 
# mod <- glm(cbind(suc, fai) ~ Growth.form + Coral.Family, miz, family=binomial)
# drop1(mod, test="Chisq")


# unique(miz$Coral.Family)
# unique(dat$Coral.Family)


# export
write.csv(dat, "data/data.csv")

