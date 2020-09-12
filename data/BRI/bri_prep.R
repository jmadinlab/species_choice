
dat<-read.csv("data/BRI/bleaching_response_index.csv")[,c(1,2)]
names(dat)<-c("species", "BRI")


head(dat)
nrow(dat)



# export
write.csv(dat, file = "data/BRI/output.csv")
##############################################################



