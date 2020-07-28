# make stereogene script

write_stereogene<-function(name_dummy_data){
  stereogene<-paste0("./Stereogene\t",
               "cfg=", name_dummy_data, ".cfg\t",
               name_dummy_data, "_track1.bedgraph\t",
               name_dummy_data, "_track2.bedgraph")
  writeLines(stereogene, "run_stereogene.sh")
}