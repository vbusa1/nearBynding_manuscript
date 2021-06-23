setwd("/Users/veronica/Documents/GitHub/RBP_binding/hold_files/alt_programs/")

# STAU2
stem<-read.delim("STAU2_HepG2/stem.out")
exterior<-read.delim("STAU2_HepG2/exterior.out")
internal<-read.delim("STAU2_HepG2/internal.out")
bulge<-read.delim("STAU2_HepG2/bulge.out")
multibranch<-read.delim("STAU2_HepG2/multibranch.out")
hairpin<-read.delim("STAU2_HepG2/hairpin.out")

png("STAU2_HepG2/visual.png",
    width = 5, height = 4, units = "in", res = 200)
plot(stem$reldist, stem$fraction, type = 'l', col = "green",
     main = "STAU2 HepG2", xlab = "Relative Distance", ylab = "Fraction",
     cex.axis = 1, cex.lab = 1, cex.main = 1.2, lwd = 2,
     ylim = c(0, .27))
lines(stem$reldist, rep(0.02, length(stem$reldist)), col =  "black")
lines(stem$reldist, bulge$fraction, col = "blue", lwd = 2)
lines(stem$reldist, multibranch$fraction, col = "red", lwd = 2)
lines(stem$reldist, hairpin$fraction, col = "purple", lwd = 2)
lines(stem$reldist, internal$fraction, col = "orange", lwd = 2)
lines(stem$reldist, exterior$fraction, col = "cyan", lwd = 2)
legend("topright",
       legend = c("Bulge", "Multibranch", "Stem",
                  "Hairpin", "Internal", "Exterior", "Random"),
       col = c("blue", "red", "green", "purple",
               "orange", "cyan", "black"), lty = 1, cex = .8, lwd = 2)
dev.off()

# PUS1
# only has 35 peaks from CLIPper so not great results
stem<-read.delim("PUS1_K562/stem.out")
exterior<-read.delim("PUS1_K562/exterior.out")
internal<-read.delim("PUS1_K562/internal.out")
bulge<-read.delim("PUS1_K562/bulge.out")
multibranch<-read.delim("PUS1_K562/multibranch.out")
hairpin<-read.delim("PUS1_K562/hairpin.out")

plot(stem$reldist, stem$fraction, type = 'l', col = "green",
     main = "PUS1 K562", xlab = "Relative Distance", ylab = "Fraction",
     cex.axis = 1, cex.lab = 1, cex.main = 1.2, lwd = 2)
lines(stem$reldist, rep(0.02, length(stem$reldist)), col =  "black")
lines(stem$reldist, bulge$fraction, col = "blue", lwd = 2)
lines(stem$reldist, multibranch$fraction, col = "red", lwd = 2)
lines(stem$reldist, hairpin$fraction, col = "purple", lwd = 2)
lines(stem$reldist, internal$fraction, col = "orange", lwd = 2)
lines(stem$reldist, exterior$fraction, col = "cyan", lwd = 2)
legend("topright",
       legend = c("Bulge", "Multibranch", "Stem",
                  "Hairpin", "Internal", "Exterior", "Random"),
       col = c("blue", "red", "green", "purple",
               "orange", "cyan", "black"), lty = 1, cex = .8, lwd = 2)

# LIN28B
# HepG2

stem<-read.delim("LIN28B_HepG2/stem.out")
exterior<-read.delim("LIN28B_HepG2/exterior.out")
internal<-read.delim("LIN28B_HepG2/internal.out")
bulge<-read.delim("LIN28B_HepG2/bulge.out")
multibranch<-read.delim("LIN28B_HepG2/multibranch.out")
hairpin<-read.delim("LIN28B_HepG2/hairpin.out")

png("LIN28B_HepG2/visual.png",
    width = 5, height = 4, units = "in", res = 200)
plot(stem$reldist, stem$fraction, type = 'l', col = "green",
     main = "LIN28B HepG2", xlab = "Relative Distance", ylab = "Fraction",
     cex.axis = 1, cex.lab = 1, cex.main = 1.2, lwd = 2,
     ylim = c(0, .27))
lines(stem$reldist, rep(0.02, length(stem$reldist)), col =  "black")
lines(stem$reldist, bulge$fraction, col = "blue", lwd = 2)
lines(stem$reldist, multibranch$fraction, col = "red", lwd = 2)
lines(stem$reldist, hairpin$fraction, col = "purple", lwd = 2)
lines(stem$reldist, internal$fraction, col = "orange", lwd = 2)
lines(stem$reldist, exterior$fraction, col = "cyan", lwd = 2)
legend("topright",
       legend = c("Bulge", "Multibranch", "Stem",
                  "Hairpin", "Internal", "Exterior", "Random"),
       col = c("blue", "red", "green", "purple",
               "orange", "cyan", "black"), lty = 1, cex = .8, lwd = 2)
dev.off()
# K562
stem<-read.delim("LIN28B_K562/stem.out")
exterior<-read.delim("LIN28B_K562/exterior.out")
internal<-read.delim("LIN28B_K562/internal.out")
bulge<-read.delim("LIN28B_K562/bulge.out")
multibranch<-read.delim("LIN28B_K562/multibranch.out")
hairpin<-read.delim("LIN28B_K562/hairpin.out")

png("LIN28B_K562/visual.png",
    width = 5, height = 4, units = "in", res = 200)
plot(stem$reldist, stem$fraction, type = 'l', col = "green",
     main = "LIN28B K562", xlab = "Relative Distance", ylab = "Fraction",
     cex.axis = 1, cex.lab = 1, cex.main = 1.2, lwd = 2,
     ylim = c(0, .27))
lines(stem$reldist, rep(0.02, length(stem$reldist)), col =  "black")
lines(stem$reldist, bulge$fraction, col = "blue", lwd = 2)
lines(stem$reldist, multibranch$fraction, col = "red", lwd = 2)
lines(stem$reldist, hairpin$fraction, col = "purple", lwd = 2)
lines(stem$reldist, internal$fraction, col = "orange", lwd = 2)
lines(stem$reldist, exterior$fraction, col = "cyan", lwd = 2)
legend("topright",
       legend = c("Bulge", "Multibranch", "Stem",
                  "Hairpin", "Internal", "Exterior", "Random"),
       col = c("blue", "red", "green", "purple",
               "orange", "cyan", "black"), lty = 1, cex = .8, lwd = 2)
dev.off()

# visualize relative positions of intronic proteins
# UKE = PRPF8
# DFO = U2AF2
# ILV = BUD13
# ZMJ = RBM22

PRPF8_RBM22<-read.delim("intronic_proteins/UKE_ZMJ.out")
U2AF2_BUD13<-read.delim("intronic_proteins/DFO_ILV.out")
U2AF2_RBM22<-read.delim("intronic_proteins/DFO_ZMJ.out")
U2AF2_PRPF8<-read.delim("intronic_proteins/DFO_UKE.out")

png("intronic_proteins/visual_PRPF8_RBM22.png",
    width = 5, height = 4, units = "in", res = 200)
plot(PRPF8_RBM22$reldist, PRPF8_RBM22$fraction, type = 'l', col = "blue",
     main = "PRPF8 vs RBM22", xlab = "Relative Distance", ylab = "Fraction",
     cex.axis = 1, cex.lab = 1, cex.main = 1.2, lwd = 2)
lines(PRPF8_RBM22$reldist, rep(0.02, length(PRPF8_RBM22$reldist)))
legend("topright",
       legend = c("random", "RBM22"),
       col = c("black", "blue"), lty = 1, cex = .8, lwd = 2)
dev.off()

png("intronic_proteins/visual_U2AF2_vs.png",
    width = 5, height = 4, units = "in", res = 200)
plot(U2AF2_BUD13$reldist, U2AF2_BUD13$fraction, type = 'l', col = "darkgreen",
     main = "U2AF2 vs BUD13, RBM22, and PRPF8", xlab = "Relative Distance", 
     ylab = "Fraction",
     cex.axis = 1, cex.lab = 1, cex.main = 1.2, lwd = 2)
lines(U2AF2_BUD13$reldist, rep(0.02, length(U2AF2_BUD13$reldist)))
lines(U2AF2_BUD13$reldist, U2AF2_RBM22$fraction, col = "blue", lwd = 2)
lines(U2AF2_PRPF8$reldist, U2AF2_PRPF8$fraction, col = "red", lwd = 2)
legend("topright",
       legend = c("random", "BUD13", "RBM22", "PRPF8"),
       col = c("black", "darkgreen", "blue", "red"), lty = 1, cex = .8, lwd = 2)
dev.off()
