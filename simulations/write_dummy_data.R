# write track 1 and track 2 bedgraphs and chromosome size file for dummy data

write_dummy_data<-function(track1, track2, name_dummy_data){
  # check for errors
  if(length(track1) != length(track2)){
    print("STOP, things aren't the same length!")
  } else{
    print("It all looks fine")
  }
  
  #write bedgraph file
  header<-paste0("#bedGraph section dummy_chr:1-",length(track1))
  track1_bedgraph<-data.frame(chr="dummy_chr",
                           start=0:(length(track1)-1),
                           stop=1:length(track1),
                           score=track1)
  track2_bedgraph<-data.frame(chr="dummy_chr",
                             start=0:(length(track1)-1),
                             stop=1:length(track1),
                             score=track2)
  writeLines(header, paste0(name_dummy_data, "_track1.bedgraph"))
  write.table(track1_bedgraph, paste0(name_dummy_data, "_track1.bedgraph"), 
              append = T,  sep="\t", row.names = F, col.names = F, quote = F)
  writeLines(header, paste0(name_dummy_data, "_track2.bedgraph"))
  write.table(track2_bedgraph, paste0(name_dummy_data, "_track2.bedgraph"), 
              append = T,  sep="\t", row.names = F, col.names = F, quote = F)
  
  # write chromosome size file
  write(paste0("dummy_chr\t", length(track1)), paste0(name_dummy_data, ".size"))
}
