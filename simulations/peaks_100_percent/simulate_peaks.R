# test what correlation looks like if there is 100% peak binding
# makes peak 100% of time

source("../write_dummy_data.R")
source("../make_config.R")
source("../write_stereogene.R")

### background peaks are double height of foreground peaks
name_dummy_data<-"bg_2_height_of_fg"
track1<-c()
track2<-c()
track3<-c()
j<-1
while(j < 50000){
  dist<-runif(1, min=5, max=100)
  height_background<-runif(1, min=4, max=200)
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

i<-1
while(i < 10000){
  dist<-runif(1, min=31, max=500)
  height_RNA<-runif(1, min=.02, max=1)
  height_protein<-height_RNA
  track2<-c(track2,
            height_RNA*c(.1,.3,1,.3,.1), #make variable-height 5 unit peak
            rep(0, dist))
  track3<-c(track3,
            height_protein*(sin(1:30/10)), #30 nt peak
            rep(0, length(track2)-length(track3)))
  i<-i+1
}
track2<-c(rep(0, length(track3)-length(track2)), track2)

track2<-track2[1:length(track1)]
track3<-track3[1:length(track1)]
track3<-track3+track1

write_dummy_data(track2, track3, name_dummy_data)
make_config(name_dummy_data)
write_stereogene(name_dummy_data)
system("sh run_stereogene.sh")

stereogene<-paste0("Stereogene\t",
                   "cfg=", name_dummy_data, ".cfg\t",
                   "background.bedgraph\t",
                   name_dummy_data, "_track2.bedgraph")
writeLines(stereogene, "run_stereogene.sh")
system("sh run_stereogene.sh")

### background peaks are five times height of foreground peaks
name_dummy_data<-"bg_5_height_of_fg"
track2<-c()
track3<-c()

i<-1
while(i < 10000){
  dist<-runif(1, min=31, max=500)
  height_RNA<-runif(1, min=.02, max=1)
  height_protein<-height_RNA
  track2<-c(track2,
            height_RNA*c(.1,.3,1,.3,.1), #make variable-height 5 unit peak
            rep(0, dist))
  track3<-c(track3,
            height_protein*(sin(1:30/10)), #30 nt peak
            rep(0, length(track2)-length(track3)))
  i<-i+1
}
track2<-c(rep(0, length(track3)-length(track2)), track2)

track2<-track2[1:length(track1)]
track3<-track3[1:length(track1)]
track3<-track3+track1

write_dummy_data(track2, track3, name_dummy_data)
make_config(name_dummy_data)
write_stereogene(name_dummy_data)
system("sh run_stereogene.sh")

stereogene<-paste0("Stereogene\t",
                   "cfg=", name_dummy_data, ".cfg\t",
                   "background.bedgraph\t",
                   name_dummy_data, "_track2.bedgraph")
writeLines(stereogene, "run_stereogene.sh")
system("sh run_stereogene.sh")

### background peaks are 20 times height of foreground peaks
name_dummy_data<-"bg_20_height_of_fg"
track2<-c()
track3<-c()

i<-1
while(i < 10000){
  dist<-runif(1, min=31, max=500)
  height_RNA<-runif(1, min=.02, max=1)
  height_protein<-height_RNA
  track2<-c(track2,
            height_RNA*c(.1,.3,1,.3,.1), #make variable-height 5 unit peak
            rep(0, dist))
  track3<-c(track3,
            height_protein*(sin(1:30/10)), #30 nt peak
            rep(0, length(track2)-length(track3)))
  i<-i+1
}
track2<-c(rep(0, length(track3)-length(track2)), track2)

track2<-track2[1:length(track1)]
track3<-track3[1:length(track1)]
track3<-track3+track1

write_dummy_data(track2, track3, name_dummy_data)
make_config(name_dummy_data)
write_stereogene(name_dummy_data)
system("sh run_stereogene.sh")

stereogene<-paste0("Stereogene\t",
                   "cfg=", name_dummy_data, ".cfg\t",
                   "background.bedgraph\t",
                   name_dummy_data, "_track2.bedgraph")
writeLines(stereogene, "run_stereogene.sh")
system("sh run_stereogene.sh")

### background peaks are same range of sizes as foreground peaks
name_dummy_data<-"bg_same_as_fg"
track2<-c()
track3<-c()

i<-1
while(i < 10000){
  dist<-runif(1, min=31, max=500)
  height_RNA<-runif(1, min=.02, max=1)
  height_protein<-height_RNA
  track2<-c(track2,
            height_RNA*c(.1,.3,1,.3,.1), #make variable-height 5 unit peak
            rep(0, dist))
  track3<-c(track3,
            height_protein*(sin(1:30/10)), #30 nt peak
            rep(0, length(track2)-length(track3)))
  i<-i+1
}
track2<-c(rep(0, length(track3)-length(track2)), track2)

track2<-track2[1:length(track1)]
track3<-track3[1:length(track1)]
track3<-track3+track1

write_dummy_data(track2, track3, name_dummy_data)
make_config(name_dummy_data)
write_stereogene(name_dummy_data)
system("sh run_stereogene.sh")

stereogene<-paste0("Stereogene\t",
                   "cfg=", name_dummy_data, ".cfg\t",
                   "background.bedgraph\t",
                   name_dummy_data, "_track2.bedgraph")
writeLines(stereogene, "run_stereogene.sh")
system("sh run_stereogene.sh")

### background peaks are 1/5 height of foreground peaks
name_dummy_data<-"bg_.2_height_of_fg"
track2<-c()
track3<-c()

i<-1
while(i < 10000){
  dist<-runif(1, min=31, max=500)
  height_RNA<-runif(1, min=.02, max=1)
  height_protein<-height_RNA
  track2<-c(track2,
            height_RNA*c(.1,.3,1,.3,.1), #make variable-height 5 unit peak
            rep(0, dist))
  track3<-c(track3,
            height_protein*(sin(1:30/10)), #30 nt peak
            rep(0, length(track2)-length(track3)))
  i<-i+1
}
track2<-c(rep(0, length(track3)-length(track2)), track2)

track2<-track2[1:length(track1)]
track3<-track3[1:length(track1)]
track3<-track3+track1

write_dummy_data(track2, track3, name_dummy_data)
make_config(name_dummy_data)
write_stereogene(name_dummy_data)
system("sh run_stereogene.sh")

stereogene<-paste0("Stereogene\t",
                   "cfg=", name_dummy_data, ".cfg\t",
                   "background.bedgraph\t",
                   name_dummy_data, "_track2.bedgraph")
writeLines(stereogene, "run_stereogene.sh")
system("sh run_stereogene.sh")

### background peaks are 1/100 height of foreground peaks
name_dummy_data<-"bg_.01_height_of_fg"
track2<-c()
track3<-c()

i<-1
while(i < 10000){
  dist<-runif(1, min=31, max=500)
  height_RNA<-runif(1, min=.02, max=1)
  height_protein<-height_RNA
  track2<-c(track2,
            height_RNA*c(.1,.3,1,.3,.1), #make variable-height 5 unit peak
            rep(0, dist))
  track3<-c(track3,
            height_protein*(sin(1:30/10)), #30 nt peak
            rep(0, length(track2)-length(track3)))
  i<-i+1
}
track2<-c(rep(0, length(track3)-length(track2)), track2)

track2<-track2[1:length(track1)]
track3<-track3[1:length(track1)]
track3<-track3+track1

write_dummy_data(track2, track3, name_dummy_data)
make_config(name_dummy_data)
write_stereogene(name_dummy_data)
system("sh run_stereogene.sh")

stereogene<-paste0("Stereogene\t",
                   "cfg=", name_dummy_data, ".cfg\t",
                  "background.bedgraph\t",
                   name_dummy_data, "_track2.bedgraph")
writeLines(stereogene, "run_stereogene.sh")
system("sh run_stereogene.sh")

### background peaks are 1/20 height of foreground peaks
name_dummy_data<-"bg_.05_height_of_fg"
track2<-c()
track3<-c()

i<-1
while(i < 10000){
  dist<-runif(1, min=31, max=500)
  height_RNA<-runif(1, min=.02, max=1)
  height_protein<-height_RNA
  track2<-c(track2,
            height_RNA*c(.1,.3,1,.3,.1), #make variable-height 5 unit peak
            rep(0, dist))
  track3<-c(track3,
            height_protein*(sin(1:30/10)), #30 nt peak
            rep(0, length(track2)-length(track3)))
  i<-i+1
}
track2<-c(rep(0, length(track3)-length(track2)), track2)

track2<-track2[1:length(track1)]
track3<-track3[1:length(track1)]
track3<-track3+track1

write_dummy_data(track2, track3, name_dummy_data)
make_config(name_dummy_data)
write_stereogene(name_dummy_data)
system("sh run_stereogene.sh")

stereogene<-paste0("Stereogene\t",
                   "cfg=", name_dummy_data, ".cfg\t",
                   "background.bedgraph\t",
                   name_dummy_data, "_track2.bedgraph")
writeLines(stereogene, "run_stereogene.sh")
system("sh run_stereogene.sh")

### background peaks are 1/40 height of foreground peaks
name_dummy_data<-"bg_.025_height_of_fg"
track2<-c()
track3<-c()
j<-1

i<-1
while(i < 10000){
  dist<-runif(1, min=31, max=500)
  height_RNA<-runif(1, min=.02, max=1)
  height_protein<-height_RNA
  track2<-c(track2,
            height_RNA*c(.1,.3,1,.3,.1), #make variable-height 5 unit peak
            rep(0, dist))
  track3<-c(track3,
            height_protein*(sin(1:30/10)), #30 nt peak
            rep(0, length(track2)-length(track3)))
  i<-i+1
}
track2<-c(rep(0, length(track3)-length(track2)), track2)

track2<-track2[1:length(track1)]
track3<-track3[1:length(track1)]
track3<-track3+track1

write_dummy_data(track2, track3, name_dummy_data)
make_config(name_dummy_data)
write_stereogene(name_dummy_data)
system("sh run_stereogene.sh")

stereogene<-paste0("Stereogene\t",
                   "cfg=", name_dummy_data, ".cfg\t",
                   "background.bedgraph\t",
                   name_dummy_data, "_track2.bedgraph")
writeLines(stereogene, "run_stereogene.sh")
system("sh run_stereogene.sh")

#### plot all together
same<-read.table("bg_same_as_fg_track1~bg_same_as_fg_track2.dist", header = T)
hundred<-read.table("bg_.01_height_of_fg_track1~bg_.01_height_of_fg_track2.dist", header = T)
five<-read.table("bg_.2_height_of_fg_track1~bg_.2_height_of_fg_track2.dist", header = T)
forty<-read.table("bg_.025_height_of_fg_track1~bg_.025_height_of_fg_track2.dist", header = T)
twenty<-read.table("bg_.05_height_of_fg_track1~bg_.05_height_of_fg_track2.dist", header = T)
half<-read.table("bg_2_height_of_fg_track1~bg_2_height_of_fg_track2.dist", header = T)
fifth<-read.table("bg_5_height_of_fg_track1~bg_5_height_of_fg_track2.dist", header = T)
twentieth<-read.table("bg_20_height_of_fg_track1~bg_20_height_of_fg_track2.dist", header = T)

pdf("all_peaks.pdf", height = 3, width = 4)
#  create the plot 
old.par <- par( no.readonly = TRUE ) 
par(oma = c( 0, 0, 0, 0 ),mar=c(3,3,3,1),mgp=c(1.6,0.45,0))
plot(hundred$x, hundred$Fg, type = 'l',col = 'black', xlim = c(-10,40),
     xlab = 'Distance',
     ylab ='Density*100',cex.axis = 0.8,  cex.lab = .9, lwd=2)
abline(v=0, col="grey", lty=2)
lines(same$x, same$Fg , col='blue',lwd=2)
lines(forty$x, forty$Fg , col='red',lwd=2)
lines(twenty$x, twenty$Fg , col='cyan',lwd=2)
lines(five$x, five$Fg , col='green',lwd=2)
lines(half$x, half$Fg , col='orange',lwd=2)
lines(fifth$x, fifth$Fg , col='purple',lwd=2)
lines(twentieth$x, twentieth$Fg , col='forestgreen',lwd=2)
legend("right", legend=c("100", "40", "20", "5", "1", "0.5", "0.2", "0.05"),
       col=c("black", "red", "cyan", "green", "blue", "orange", "purple", "forestgreen"), lty=1, cex=0.7,
       title = "foreground/background\nsignal ratio", bty = "n")
dev.off()

marks<-c(0, .1, 1, 10, 100)
pdf("all_points.pdf", height = 4, width = 4)
#  create the plot 
old.par <- par( no.readonly = TRUE ) 
par(oma = c( 0, 0, 0, 0 ),mar=c(3,3,3,1),mgp=c(1.6,0.45,0))
plot(c(1, 5, 20, 40, 100, .5, .2, .05),
     log = "x",
     xaxt="n",
     c(max(same$Fg),
       max(five$Fg),
       max(twenty$Fg),
       max(forty$Fg),
       max(hundred$Fg),
       max(half$Fg),
       max(fifth$Fg),
       max(twentieth$Fg)),
     xlab="Ratio of Forground/Background",
     ylab="Maximum Density * 100",
     pch=19, lwd=5)
axis(1, at=marks,labels=marks)
dev.off()
