# test what correlation looks like if peaks are uneven widths across transcriptome

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

pos<-rbeta(10000, 5, 5) * length(track1)
track2<-rep(0, length(track1))
track3<-rep(0, length(track1))
track4<-rep(0, length(track1))
track5<-rep(0, length(track1))
track6<-rep(0, length(track1))
track7<-rep(0, length(track1))
for(n in pos){
  height_protein<-height_RNA<-runif(1, min=.02, max=1)
  track2[n:(n+4)]<-track2[n:(n+4)]+(height_RNA*c(.1,.3,1,.3,.1)) #make variable-height 5 unit peak
  track3[n:(n+29)]<-track3[n:(n+29)]+(height_protein*(sin(1:30/(30/pi)))) #make variable-height 30 unit peak
  width_protein_4<-round(runif(1, min=28, max=32), 0)
  if(width_protein_4 %% 2 == 0){
    i <- width_protein_4/2
    track4[(n+(15-i)):(n+(14+i))]<-track4[(n+(15-i)):(n+(14+i))]+
      (height_protein*(sin(1:width_protein_4/(width_protein_4/pi))))
  }else{
    i <- width_protein_4/2 + .5
    track4[(n+(15-i)):(n+(13+i))]<-track4[(n+(15-i)):(n+(13+i))]+
      (height_protein*(sin(1:width_protein_4/(width_protein_4/pi))))
  }
  width_protein_5<-round(runif(1, min=25, max=35), 0)
  if(width_protein_5 %% 2 == 0){
    i <- width_protein_5/2
    track5[(n+(15-i)):(n+(14+i))]<-track5[(n+(15-i)):(n+(14+i))]+
      (height_protein*(sin(1:width_protein_5/(width_protein_5/pi))))
  }else{
    i <- width_protein_5/2 + .5
    track5[(n+(15-i)):(n+(13+i))]<-track5[(n+(15-i)):(n+(13+i))]+
      (height_protein*(sin(1:width_protein_5/(width_protein_5/pi))))
  }
  width_protein_6<-round(runif(1, min=20, max=40), 0)
  if(width_protein_6 %% 2 == 0){
    i <- width_protein_6/2
    track6[(n+(15-i)):(n+(14+i))]<-track6[(n+(15-i)):(n+(14+i))]+
      (height_protein*(sin(1:width_protein_6/(width_protein_6/pi))))
  }else{
    i <- width_protein_6/2 + .5
    track6[(n+(15-i)):(n+(13+i))]<-track6[(n+(15-i)):(n+(13+i))]+
      (height_protein*(sin(1:width_protein_6/(width_protein_6/pi))))
  }
  width_protein_7<-round(runif(1, min=10, max=50), 0)
  if(width_protein_7 %% 2 == 0){
    i <- width_protein_7/2
    track7[(n+(15-i)):(n+(14+i))]<-track7[(n+(15-i)):(n+(14+i))]+
      (height_protein*(sin(1:width_protein_7/(width_protein_7/pi))))
  }else{
    i <- width_protein_7/2 + .5
    track7[(n+(15-i)):(n+(13+i))]<-track7[(n+(15-i)):(n+(13+i))]+
      (height_protein*(sin(1:width_protein_7/(width_protein_7/pi))))
  }
}
track3<-track3+track1
track4<-track4+track1
track5<-track5+track1
track6<-track6+track1
track7<-track7+track1

write_dummy_data(track2, track3, "peak_30")
make_config("peak_30")
write_stereogene("peak_30")
system("sh run_stereogene.sh")

stereogene<-paste0("Stereogene\t",
                   "cfg=peak_30.cfg\t", # all tracks will have same size so doesn't matter which is used
                   "peak_30_track1.bedgraph\t", # RNA track is constant across all comparisons
                   "background.bedgraph")
writeLines(stereogene, "run_stereogene.sh")
system("sh run_stereogene.sh")

write_dummy_data(track2, track4, "peak_28_32")
make_config("peak_28_32")
write_stereogene("peak_28_32")
system("sh run_stereogene.sh")

write_dummy_data(track2, track5, "peak_25_35")
make_config("peak_25_35")
write_stereogene("peak_25_35")
system("sh run_stereogene.sh")

write_dummy_data(track2, track6, "peak_20_40")
make_config("peak_20_40")
write_stereogene("peak_20_40")
system("sh run_stereogene.sh")

write_dummy_data(track2, track7, "peak_10_50")
make_config("peak_10_50")
write_stereogene("peak_10_50")
system("sh run_stereogene.sh")

#### plot all together
background<-read.table("peak_30_track1~background.dist", header = T)
peaks_30<-read.table("peak_30_track1~peak_30_track2.dist", header = T)
peaks_28_32<-read.table("peak_28_32_track1~peak_28_32_track2.dist", header = T)
peaks_25_35<-read.table("peak_25_35_track1~peak_25_35_track2.dist", header = T)
peaks_20_40<-read.table("peak_20_40_track1~peak_20_40_track2.dist", header = T)
peaks_10_50<-read.table("peak_10_50_track1~peak_10_50_track2.dist", header = T)

peaks_30$Fg<-peaks_30$Fg - background$Fg
peaks_30$Bkg<-peaks_30$Bkg - background$Bkg
peaks_28_32$Fg<-peaks_28_32$Fg - background$Fg
peaks_28_32$Bkg<-peaks_28_32$Bkg - background$Bkg
peaks_25_35$Fg<-peaks_25_35$Fg - background$Fg
peaks_25_35$Bkg<-peaks_25_35$Bkg - background$Bkg
peaks_20_40$Fg<-peaks_20_40$Fg - background$Fg
peaks_20_40$Bkg<-peaks_20_40$Bkg - background$Bkg
peaks_10_50$Fg<-peaks_10_50$Fg - background$Fg
peaks_10_50$Bkg<-peaks_10_50$Bkg - background$Bkg

pdf("compare_peaks_widths.pdf", height = 3.7, width = 3.7)
#  create the plot
old.par <- par( no.readonly = TRUE )
par(oma = c( 0, 0, 0, 0 ),
    mar=c(3,3,3,1),
    mgp=c(1.6,0.45,0))
plot(-5000:4999,
     peaks_30$Bkg,
     type = 'l', col = 'black',
     xlim = c(-10,40),
     ylim = c(0,18),
     xlab = 'Distance',
     ylab ='Density*100',
     cex.axis = 0.8,
     cex.lab = .9,
     lwd=2)
lines(-5000:4999,
      peaks_28_32$Bkg,
      col = "black", lwd = 2)
lines(-5000:4999,
      peaks_25_35$Bkg,
      col = "black", lwd = 2)
lines(-5000:4999,
      peaks_20_40$Bkg,
      col = "black", lwd = 2)
lines(-5000:4999,
      peaks_10_50$Bkg,
      col = "black", lwd = 2)
abline(v=0, col="grey", lty=2)
lines(-5000:4999,
      peaks_30$Fg,
      col='blue',lwd=2)
lines(-5000:4999,
      peaks_28_32$Fg,
      col='red',lwd=2)
lines(-5000:4999,
      peaks_25_35$Fg,
      col='forestgreen',lwd=2)
lines(-5000:4999,
      peaks_20_40$Fg,
      col='goldenrod1',lwd=2)
lines(-5000:4999,
      peaks_10_50$Fg,
      col='cadetblue',lwd=2)
legend("right", legend=c("background", "30", "28-32", "25-35", "20-40", "10-50"),
       col=c("black", "blue", "red", "forestgreen", "goldenrod1", "cadetblue"), lty=1, cex=0.7,
       title = "Peak Width Range", bty = "n")
dev.off()
