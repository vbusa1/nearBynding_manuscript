# test what correlation looks like if peaks are unevenly distributed across transcriptome

source("write_dummy_data.R")
source("make_config.R")
source("write_stereogene.R")

track1<-c() #background track
j<-1
while(j < 50000){
  dist<-runif(1, min=5, max=100)
  height_background<-runif(1, min=.02, max=1)
  track1<-c(track1,
            height_background*(sin(1:30/10)),
            rep(0, dist))
  j<-j+1
}

header<-paste0("#bedGraph section dummy_chr:1-",length(track1))
track1_bedgraph<-data.frame(chr="dummy_chr",
                            start=0:(length(track1)-1),
                            stop=1:length(track1),
                            score=track1)
writeLines(header, "background.bedgraph")
write.table(track1_bedgraph, "background.bedgraph", append = T,  sep="\t", row.names = F, col.names = F, quote = F)


name_dummy_data<-"even_dist"
pos<-rbeta(10000, 5, 5) * length(track1)
track2<-rep(0, length(track1))
track3<-rep(0, length(track1))
for(n in pos){
  height_protein<-height_RNA<-runif(1, min=.02, max=1)
  track2[n:(n+4)]<-track2[n:(n+4)]+(height_RNA*c(.1,.3,1,.3,.1)) #make variable-height 5 unit peak
  track3[n:(n+29)]<-track3[n:(n+29)]+(height_protein*(sin(1:30/10))) #make variable-height 30 unit peak
}
track3<-track3[1:length(track1)]+track1

write_dummy_data(track2, track3, name_dummy_data)
make_config(name_dummy_data)
write_stereogene(name_dummy_data)
system("sh run_stereogene.sh")

stereogene<-paste0("Stereogene\t",
                   "cfg=", name_dummy_data, ".cfg\t",
                   name_dummy_data, "_track1.bedgraph\t",
                   "background.bedgraph")
writeLines(stereogene, "run_stereogene.sh")
system("sh run_stereogene.sh")

name_dummy_data<-"skew_left_1"
pos<-rbeta(10000, 1, 5) * length(track1)
track2<-rep(0, length(track1))
track3<-rep(0, length(track1))
for(n in pos){
  height_protein<-height_RNA<-runif(1, min=.02, max=1)
  track2[n:(n+4)]<-track2[n:(n+4)]+(height_RNA*c(.1,.3,1,.3,.1)) #make variable-height 5 unit peak
  track3[n:(n+29)]<-track3[n:(n+29)]+(height_protein*(sin(1:30/10))) #make variable-height 30 unit peak
}
track3<-track3[1:length(track1)]+track1

write_dummy_data(track2, track3, name_dummy_data)
make_config(name_dummy_data)
write_stereogene(name_dummy_data)
system("sh run_stereogene.sh")

stereogene<-paste0("Stereogene\t",
                   "cfg=", name_dummy_data, ".cfg\t",
                   name_dummy_data, "_track1.bedgraph\t",
                   "background.bedgraph")
writeLines(stereogene, "run_stereogene.sh")
system("sh run_stereogene.sh")


name_dummy_data<-"skew_left_3"
pos<-rbeta(10000, 3, 5) * length(track1)
track2<-rep(0, length(track1))
track3<-rep(0, length(track1))
for(n in pos){
  height_protein<-height_RNA<-runif(1, min=.02, max=1)
  track2[n:(n+4)]<-track2[n:(n+4)]+(height_RNA*c(.1,.3,1,.3,.1)) #make variable-height 5 unit peak
  track3[n:(n+29)]<-track3[n:(n+29)]+(height_protein*(sin(1:30/10))) #make variable-height 30 unit peak
}
track3<-track3[1:length(track1)]+track1

write_dummy_data(track2, track3, name_dummy_data)
make_config(name_dummy_data)
write_stereogene(name_dummy_data)
system("sh run_stereogene.sh")

stereogene<-paste0("Stereogene\t",
                   "cfg=", name_dummy_data, ".cfg\t",
                   name_dummy_data, "_track1.bedgraph\t",
                   "background.bedgraph")
writeLines(stereogene, "run_stereogene.sh")
system("sh run_stereogene.sh")


name_dummy_data<-"skew_right_1"
pos<-rbeta(10000, 5, 1) * length(track1)
track2<-rep(0, length(track1))
track3<-rep(0, length(track1))
for(n in pos){
  height_protein<-height_RNA<-runif(1, min=.02, max=1)
  track2[n:(n+4)]<-track2[n:(n+4)]+(height_RNA*c(.1,.3,1,.3,.1)) #make variable-height 5 unit peak
  track3[n:(n+29)]<-track3[n:(n+29)]+(height_protein*(sin(1:30/10))) #make variable-height 30 unit peak
}
track3<-track3[1:length(track1)]+track1

write_dummy_data(track2, track3, name_dummy_data)
make_config(name_dummy_data)
write_stereogene(name_dummy_data)
system("sh run_stereogene.sh")

stereogene<-paste0("Stereogene\t",
                   "cfg=", name_dummy_data, ".cfg\t",
                   name_dummy_data, "_track1.bedgraph\t",
                   "background.bedgraph")
writeLines(stereogene, "run_stereogene.sh")
system("sh run_stereogene.sh")


name_dummy_data<-"skew_right_3"
pos<-rbeta(10000, 5, 3) * length(track1)
track2<-rep(0, length(track1))
track3<-rep(0, length(track1))
for(n in pos){
  height_protein<-height_RNA<-runif(1, min=.02, max=1)
  track2[n:(n+4)]<-track2[n:(n+4)]+(height_RNA*c(.1,.3,1,.3,.1)) #make variable-height 5 unit peak
  track3[n:(n+29)]<-track3[n:(n+29)]+(height_protein*(sin(1:30/10))) #make variable-height 30 unit peak
}
track3<-track3[1:length(track1)]+track1

write_dummy_data(track2, track3, name_dummy_data)
make_config(name_dummy_data)
write_stereogene(name_dummy_data)
system("sh run_stereogene.sh")

stereogene<-paste0("Stereogene\t",
                   "cfg=", name_dummy_data, ".cfg\t",
                   name_dummy_data, "_track1.bedgraph\t",
                   "background.bedgraph")
writeLines(stereogene, "run_stereogene.sh")
system("sh run_stereogene.sh")

#### plot all together
same<-read.table("even_dist_track1~even_dist_track2.dist", header = T)
left_1<-read.table("skew_left_1_track1~skew_left_1_track2.dist", header = T)
left_3<-read.table("skew_left_3_track1~skew_left_3_track2.dist", header = T)
right_1<-read.table("skew_right_1_track1~skew_right_1_track2.dist", header = T)
right_3<-read.table("skew_right_3_track1~skew_right_3_track2.dist", header = T)

files<-c("even_dist", "skew_left_1", "skew_left_3", "skew_right_1", "skew_right_3")
data<-as.data.frame(matrix(NA, nrow = 10, ncol = 10002))
colnames(data)<-c("skew", "corr", -5000:4999)
n<-1
for(skew in files){
    data[c(n, n+1),c("skew", "corr")]<-c(skew, skew, "Fg", "Bkg")
    data[n, 3:10002]<- read.table(paste0(skew, "_track1~", skew, "_track2.dist"),header = T)[,3]
    data[n+1, 3:10002]<- read.table(paste0(skew, "_track1~", "background.dist"),header = T)[,3]
    n<-n+2
}

pdf("compare_skew_peaks.pdf", height = 3.7, width = 3.7)
#  create the plot
old.par <- par( no.readonly = TRUE )
par(oma = c( 0, 0, 0, 0 ),
    mar=c(3,3,3,1),
    mgp=c(1.6,0.45,0))
plot(-5000:4999, 
     (data[which(data$skew == "even_dist"  & data$corr == "Fg"), 3:10002]-
        data[which(data$skew == "even_dist"  & data$corr == "Bkg"), 3:10002]),
     type = 'l', col = 'black',
     xlim = c(-10,40),
     xlab = 'Distance',
     ylab ='Density*100',
     cex.axis = 0.8,
     cex.lab = .9,
     lwd=2)
abline(v=0, col="grey", lty=2)
lines(-5000:4999, 
      (data[which(data$skew == "skew_left_1" & data$corr == "Fg"), 3:10002]-
         data[which(data$skew == "skew_left_1" & data$corr == "Bkg"), 3:10002]),
      col='blue',lwd=2)
lines(-5000:4999, 
      (data[which(data$skew == "skew_left_3" & data$corr == "Fg"), 3:10002]-
         data[which(data$skew == "skew_left_3" & data$corr == "Bkg"), 3:10002]),
      col='red',lwd=2)
lines(-5000:4999, 
      (data[which(data$skew == "skew_right_1" & data$corr == "Fg"), 3:10002]-
         data[which(data$skew == "skew_right_1" & data$corr == "Bkg"), 3:10002]),
      col='goldenrod1',lwd=2)
lines(-5000:4999, 
      (data[which(data$skew == "skew_right_3" & data$corr == "Fg"), 3:10002]-
         data[which(data$skew == "skew_right_3" & data$corr == "Bkg"), 3:10002]),
      col='forestgreen',lwd=2)
legend("right", legend=c("None", "Far Left", "Left", "Far Right", "Right"),
       col=c("black", "blue", "red", "goldenrod1", "forestgreen"), lty=1, cex=0.7,
       title = "Peak Distribution Skew", bty = "n")
dev.off()
