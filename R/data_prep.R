

# Coral Traits Database means  - www.coraltraits.org
#source("data/traitdatabase/ctb_prep.R")
ctb<-read.csv("data/traitdatabase/output.csv")
ctb<-ctb[ctb$Zooxanthellate=="zooxanthellate",] # zooxanthellate only
ctb<-ctb[!is.na(ctb$Abundance.GBR),] # gbr only
nrow(ctb)


# trait biogeography - McWilliam et al. 2018 PNAS
dat<-read.csv("data/pnas2018/traitbiogeography.csv")
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


# Families
tax<-read.csv("data/taxonomy.csv")
#head(tax)








# export
write.csv(dat, "data/data.csv")

