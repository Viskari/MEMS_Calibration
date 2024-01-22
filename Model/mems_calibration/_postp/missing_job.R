#MISSING JOB ON 2000 JOB LUNCHED

err<-list.files('/eos/jeodpp/data/projects/SOIL-NACA/MEMS/mems_lucas/_out', pattern = "MAOM*")
err<-gsub('MAOM', '', err)
err<-gsub('.tif', '', err)
err<-data.frame("ID"= as.numeric(err), "VAL"=1)

err2<-data.frame("ID"=seq(1:2000))

missSOC<-merge(err2, err, by.x="ID", all=TRUE)

missSOC<-subset(missSOC, is.na(VAL))            


write.csv(missSOC, "/eos/jeodpp/data/projects/SOIL-NACA/MEMS/mems_lucas/_postp/miss_job.csv", row.names = F, quote = F)

print(dim(missSOC))
