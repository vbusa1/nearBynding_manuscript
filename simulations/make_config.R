# create a config file for use by StereoGene

make_config<-function(name_dummy_data, wSize=10000, wStep=10000, bin=1){
  conf<-paste("#!bash/bin", "",
              "profPath =./",
              "trackPath=./",
              "resPath=./",
              paste0("chrom=", name_dummy_data, ".size"),
              "NA=1",
              "Rscript=1",
              "outRes=TAB",
              "writeDistr=NONE",
              "Distances=0",
              paste0("bin=", bin),
              paste0("wSize=", wSize),
              paste0("wStep=", wStep),
              "maxZero=99",
              "maxNA=99",
              "nShuffle=100000",
              "outLC=1",
              sep="\n")
  writeLines(conf, paste0(name_dummy_data, ".cfg"))
}