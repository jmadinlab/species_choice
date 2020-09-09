
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




# export
write.csv(dat, file = "data/pnas2018/output.csv")
##############################################################



