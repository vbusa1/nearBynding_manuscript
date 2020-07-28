files<-c("same_as", ".01_height_of",".2_height_of",".025_height_of",".05_height_of","2_height_of","5_height_of","20_height_of")
data<-as.data.frame(matrix(NA, nrow = 64, ncol = 10003))
colnames(data)<-c("concordance", "ratio", "corr", -5000:4999)
n<-1
for(concordance in c(40, 60, 80, 100)){
        for(ratio in c(.05, .2, .5, 1, 5, 20, 40, 100)){
                if(ratio == .05){file<-files[8]}
                else if(ratio == 40){file<-files[4]}
                else if(ratio == .2){file<-files[7]}
                else if(ratio == .5){file<-files[6]}
                else if(ratio == 1){file<-files[1]}
                else if(ratio == 5){file<-files[3]}
                else if(ratio == 20){file<-files[5]}
                else if(ratio == 100){file<-files[2]}
                data[c(n, n+1),c("concordance", "ratio", "corr")]<-c(concordance, concordance,
                                                                     ratio, ratio, "Fg", "Bkg")
                data[n, 4:10003]<- read.table(paste0("peaks_", concordance, "_percent/bg_",
                                                     file, "_fg_track1~bg_", file, "_fg_track2.dist"),header = T)[,3]
                data[n+1, 4:10003]<- read.table(paste0("peaks_", concordance, "_percent/bg_",
                                                     file, "_fg_track1~background.dist"),header = T)[,3]
                n<-n+2
        }
}

pdf("compare_peaks_concordance.pdf", height = 3.7, width = 3.7)
#  create the plot
old.par <- par( no.readonly = TRUE )
par(oma = c( 0, 0, 0, 0 ),
    mar=c(3,3,3,1),
    mgp=c(1.6,0.45,0))
plot(-5000:4999, 
     (data[which(data$concordance == 100 & data$ratio == 100 & data$corr == "Fg"), 4:10003]-
              data[which(data$concordance == 100 & data$ratio == 100 & data$corr == "Bkg"), 4:10003]),
     type = 'l', col = 'black',
     xlim = c(-10,40),
     xlab = 'Distance',
     ylab ='Density x 100',
     cex.axis = 0.8,
     cex.lab = .9,
     lwd=2)
abline(v=0, col="grey", lty=2)
lines(-5000:4999, 
      (data[which(data$concordance == 80 & data$ratio == 100 & data$corr == "Fg"), 4:10003]-
               data[which(data$concordance == 80 & data$ratio == 100 & data$corr == "Bkg"), 4:10003]),
      col='blue',lwd=2)
lines(-5000:4999, 
      (data[which(data$concordance == 60 & data$ratio == 100 & data$corr == "Fg"), 4:10003]-
               data[which(data$concordance == 60 & data$ratio == 100 & data$corr == "Bkg"), 4:10003]),
      col='red',lwd=2)
lines(-5000:4999, 
      (data[which(data$concordance == 40 & data$ratio == 100 & data$corr == "Fg"), 4:10003]-
               data[which(data$concordance == 40 & data$ratio == 100 & data$corr == "Bkg"), 4:10003]),
      col='goldenrod1',lwd=2)
legend("right", legend=c("100", "80", "60", "40"),
       col=c("black", "blue", "red", "goldenrod1"), lty=1, cex=0.7,
       title = "peak concordance\nbetween tracks", bty = "n")
dev.off()

pdf("compare_peaks_background.pdf", height = 3.7, width = 3.7)
#  create the plot
old.par <- par( no.readonly = TRUE )
par(oma = c( 0, 0, 0, 0 ),
    mar=c(3,3,3,1),
    mgp=c(1.6,0.45,0))
plot(-5000:4999, 
     (data[which(data$concordance == 100 & data$ratio == 100 & data$corr == "Fg"), 4:10003]-
              data[which(data$concordance == 100 & data$ratio == 100 & data$corr == "Bkg"), 4:10003]),
     type = 'l', col = 'black',
     xlim = c(-10,40),
     xlab = 'Distance',
     ylab ='Density x 100',
     cex.axis = 0.8,
     cex.lab = .9,
     lwd=2)
abline(v=0, col="grey", lty=2)
lines(-5000:4999, 
      (data[which(data$concordance == 100 & data$ratio == 40 & data$corr == "Fg"), 4:10003]-
               data[which(data$concordance == 100 & data$ratio == 40 & data$corr == "Bkg"), 4:10003]),
      col='blue',lwd=2)
lines(-5000:4999, 
      (data[which(data$concordance == 100 & data$ratio == 20 & data$corr == "Fg"), 4:10003]-
               data[which(data$concordance == 100 & data$ratio == 20 & data$corr == "Bkg"), 4:10003]),
      col='red',lwd=2)
lines(-5000:4999, 
      (data[which(data$concordance == 100 & data$ratio == 5 & data$corr == "Fg"), 4:10003]-
               data[which(data$concordance == 100 & data$ratio == 5 & data$corr == "Bkg"), 4:10003]),
      col='goldenrod1',lwd=2)
lines(-5000:4999, 
      (data[which(data$concordance == 100 & data$ratio == 1 & data$corr == "Fg"), 4:10003]-
               data[which(data$concordance == 100 & data$ratio == 1 & data$corr == "Bkg"), 4:10003]),
      col='forestgreen',lwd=2)
legend("right", legend=c("1", "0.4", "0.2", "0.05", "0.01"),
       col=c("black", "blue", "red", "goldenrod1", "forestgreen"), lty=1, cex=0.7,
       title = "ratio of foreground\nto background", bty = "n")
dev.off()

marks<-c(0, .001, .01, .1, 1)
pdf("compare_all_points.pdf", height = 3.7, width = 3.7)
#  create the plot 
old.par <- par( no.readonly = TRUE ) 
par(oma = c( 0, 0, 0, 0 ),mar=c(3,3,3,1),mgp=c(1.6,0.45,0))
plot(c(.0005, .002, .005, .01, .05, .2, .4, 1),
     log = "x",
     xaxt="n",
     c(max(data[which(data$concordance == 100 & data$ratio == .05 & data$corr == "Fg"), 4:10003]-
                    data[which(data$concordance == 100 & data$ratio == .05 & data$corr == "Bkg"), 4:10003]),
       max(data[which(data$concordance == 100 & data$ratio == .2 & data$corr == "Fg"), 4:10003]-
                   data[which(data$concordance == 100 & data$ratio == .2 & data$corr == "Bkg"), 4:10003]),
       max(data[which(data$concordance == 100 & data$ratio == .5 & data$corr == "Fg"), 4:10003]-
                   data[which(data$concordance == 100 & data$ratio == .5 & data$corr == "Bkg"), 4:10003]),
       max(data[which(data$concordance == 100 & data$ratio == 1 & data$corr == "Fg"), 4:10003]-
                   data[which(data$concordance == 100 & data$ratio == 1 & data$corr == "Bkg"), 4:10003]),
       max(data[which(data$concordance == 100 & data$ratio == 5 & data$corr == "Fg"), 4:10003]-
                   data[which(data$concordance == 100 & data$ratio == 5 & data$corr == "Bkg"), 4:10003]),
       max(data[which(data$concordance == 100 & data$ratio == 20 & data$corr == "Fg"), 4:10003]-
                   data[which(data$concordance == 100 & data$ratio == 20 & data$corr == "Bkg"), 4:10003]),
       max(data[which(data$concordance == 100 & data$ratio == 40 & data$corr == "Fg"), 4:10003]-
                   data[which(data$concordance == 100 & data$ratio == 40 & data$corr == "Bkg"), 4:10003]),
       max(data[which(data$concordance == 100 & data$ratio == 100 & data$corr == "Fg"), 4:10003]-
                   data[which(data$concordance == 100 & data$ratio == 100 & data$corr == "Bkg"), 4:10003])),
     xlab="Ratio of Foreground/Background",
     ylab="Maximum Density * 100",
     pch=19, lwd=5)
points(c(.0005, .002, .005, .01, .05, .2, .4, 1),
       c(max(data[which(data$concordance == 80 & data$ratio == .05 & data$corr == "Fg"), 4:10003]-
                     data[which(data$concordance == 80 & data$ratio == .05 & data$corr == "Bkg"), 4:10003]),
         max(data[which(data$concordance == 80 & data$ratio == .2 & data$corr == "Fg"), 4:10003]-
                     data[which(data$concordance == 80 & data$ratio == .2 & data$corr == "Bkg"), 4:10003]),
         max(data[which(data$concordance == 80 & data$ratio == .5 & data$corr == "Fg"), 4:10003]-
                     data[which(data$concordance == 80 & data$ratio == .5 & data$corr == "Bkg"), 4:10003]),
         max(data[which(data$concordance == 80 & data$ratio == 1 & data$corr == "Fg"), 4:10003]-
                     data[which(data$concordance == 80 & data$ratio == 1 & data$corr == "Bkg"), 4:10003]),
         max(data[which(data$concordance == 80 & data$ratio == 5 & data$corr == "Fg"), 4:10003]-
                     data[which(data$concordance == 80 & data$ratio == 5 & data$corr == "Bkg"), 4:10003]),
         max(data[which(data$concordance == 80 & data$ratio == 20 & data$corr == "Fg"), 4:10003]-
                     data[which(data$concordance == 80 & data$ratio == 20 & data$corr == "Bkg"), 4:10003]),
         max(data[which(data$concordance == 80 & data$ratio == 40 & data$corr == "Fg"), 4:10003]-
                     data[which(data$concordance == 80 & data$ratio == 40 & data$corr == "Bkg"), 4:10003]),
         max(data[which(data$concordance == 80 & data$ratio == 100 & data$corr == "Fg"), 4:10003]-
                     data[which(data$concordance == 80 & data$ratio == 100 & data$corr == "Bkg"), 4:10003])),
       pch=19, lwd=5, col="red")
points(c(.0005, .002, .005, .01, .05, .2, .4, 1),
       c(max(data[which(data$concordance == 60 & data$ratio == .05 & data$corr == "Fg"), 4:10003]-
                     data[which(data$concordance == 60 & data$ratio == .05 & data$corr == "Bkg"), 4:10003]),
         max(data[which(data$concordance == 60 & data$ratio == .2 & data$corr == "Fg"), 4:10003]-
                     data[which(data$concordance == 60 & data$ratio == .2 & data$corr == "Bkg"), 4:10003]),
         max(data[which(data$concordance == 60 & data$ratio == .5 & data$corr == "Fg"), 4:10003]-
                     data[which(data$concordance == 60 & data$ratio == .5 & data$corr == "Bkg"), 4:10003]),
         max(data[which(data$concordance == 60 & data$ratio == 1 & data$corr == "Fg"), 4:10003]-
                     data[which(data$concordance == 60 & data$ratio == 1 & data$corr == "Bkg"), 4:10003]),
         max(data[which(data$concordance == 60 & data$ratio == 5 & data$corr == "Fg"), 4:10003]-
                     data[which(data$concordance == 60 & data$ratio == 5 & data$corr == "Bkg"), 4:10003]),
         max(data[which(data$concordance == 60 & data$ratio == 20 & data$corr == "Fg"), 4:10003]-
                     data[which(data$concordance == 60 & data$ratio == 20 & data$corr == "Bkg"), 4:10003]),
         max(data[which(data$concordance == 60 & data$ratio == 40 & data$corr == "Fg"), 4:10003]-
                     data[which(data$concordance == 60 & data$ratio == 40 & data$corr == "Bkg"), 4:10003]),
         max(data[which(data$concordance == 60 & data$ratio == 100 & data$corr == "Fg"), 4:10003]-
                     data[which(data$concordance == 60 & data$ratio == 100 & data$corr == "Bkg"), 4:10003])),
       pch=19, lwd=5, col="blue")
points(c(.0005, .002, .005, .01, .05, .2, .4, 1),
       c(max(data[which(data$concordance == 40 & data$ratio == .05 & data$corr == "Fg"), 4:10003]-
                     data[which(data$concordance == 40 & data$ratio == .05 & data$corr == "Bkg"), 4:10003]),
         max(data[which(data$concordance == 40 & data$ratio == .2 & data$corr == "Fg"), 4:10003]-
                     data[which(data$concordance == 40 & data$ratio == .2 & data$corr == "Bkg"), 4:10003]),
         max(data[which(data$concordance == 40 & data$ratio == .5 & data$corr == "Fg"), 4:10003]-
                     data[which(data$concordance == 40 & data$ratio == .5 & data$corr == "Bkg"), 4:10003]),
         max(data[which(data$concordance == 40 & data$ratio == 1 & data$corr == "Fg"), 4:10003]-
                     data[which(data$concordance == 40 & data$ratio == 1 & data$corr == "Bkg"), 4:10003]),
         max(data[which(data$concordance == 40 & data$ratio == 5 & data$corr == "Fg"), 4:10003]-
                     data[which(data$concordance == 40 & data$ratio == 5 & data$corr == "Bkg"), 4:10003]),
         max(data[which(data$concordance == 40 & data$ratio == 20 & data$corr == "Fg"), 4:10003]-
                     data[which(data$concordance == 40 & data$ratio == 20 & data$corr == "Bkg"), 4:10003]),
         max(data[which(data$concordance == 40 & data$ratio == 40 & data$corr == "Fg"), 4:10003]-
                     data[which(data$concordance == 40 & data$ratio == 40 & data$corr == "Bkg"), 4:10003]),
         max(data[which(data$concordance == 40 & data$ratio == 100 & data$corr == "Fg"), 4:10003]-
                     data[which(data$concordance == 40 & data$ratio == 100 & data$corr == "Bkg"), 4:10003])),
       pch=19, lwd=5, col="goldenrod1")
axis(1, at=marks,labels=marks)
legend(x=.001, y=18, legend=c("100%", "80%", "60%", "40%"),
       col=c("black", "red", "blue", "goldenrod1"), lwd=3, cex=0.8,
       title = "peak concordance\nbetween tracks", bty = "n")
dev.off()


