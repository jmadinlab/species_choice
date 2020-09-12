
library("ggplot2")
library("cowplot")

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


colnames(dat)
# TRAIT MISTAKES.. 
 #Tables should be given highest SA & space (original 2, 2, 2)
dat$cat_SA_vol[dat$raw_growth_form=="tables_or_plates"]<-5
dat$cat_spacesize[dat$raw_growth_form=="tables_or_plates"]<-4
dat$cat_colonyheight[dat$raw_growth_form=="tables_or_plates"]<-3

# branching closed? Bushy or staghorn?  (original 4, 3, 4)
#dat$cat_spacesize[dat$raw_growth_form=="branching_closed"]<-2
#dat$cat_SA_vol[dat$raw_growth_form=="branching_closed"]<-4
#dat$cat_colonyheight[dat$raw_growth_form=="branching_closed"]<-3

#dat[dat$raw_growth_form=="branching_closed","species"]

# export
write.csv(dat, file = "data/pnas2018/output.csv")
##############################################################



