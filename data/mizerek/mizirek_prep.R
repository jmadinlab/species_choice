
dat<-read.csv("data/mizerek/338_2018_1702_MOESM1_ESM.csv")[,c("Revised.species.name","BI")]
names(dat)<-c("species", "BI")

dat<-aggregate(BI~species, dat, mean)

head(dat)




# export
write.csv(dat, file = "data/mizerek/output.csv")
##############################################################



