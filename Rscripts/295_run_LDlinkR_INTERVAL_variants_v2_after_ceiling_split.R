suppressMessages(library("optparse", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("LDlinkR"))
suppressMessages(library("plyr"))
suppressMessages(library("data.table"))
suppressMessages(library("crayon"))
suppressMessages(library("withr"))
suppressMessages(library("farver"))
suppressMessages(library("dplyr"))
suppressMessages(library("withr"))



opt = NULL

options(warn = 1)

data_wrangling = function(option_list)
{

  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  path_LDlinkR = paste(out,'LDlinkR','/',sep='')
  
  if (file.exists(path_LDlinkR)){
    
    
  }else{
    
    dir.create(file.path(path_LDlinkR))
    
  }#path_LDlinkR
  
  
  #### READ covered_variants ----
  
  covered_variants<-as.data.frame(fread(file=opt$covered_variants, sep="\t", header=F), stringsAsFactors=F)
  
  colnames(covered_variants)<-'pos_for_LDlinkR'
  
  
  cat("covered_variants_0\n")
  cat(str(covered_variants))
  cat("\n")
  
  
  #### READ input_chunk ----
  
  input_chunk<-as.data.frame(fread(file=opt$input_chunk, sep="\t", header=T), stringsAsFactors=F)
  

  cat("input_chunk_0\n")
  cat(str(input_chunk))
  cat("\n")
  cat(str(unique(input_chunk$x)))
  cat("\n")
 
  
  colnames(input_chunk)[which(colnames(input_chunk) == 'x')]<-'pos_for_LDlinkR'
  
  
  cat("input_chunk_1\n")
  cat(str(input_chunk))
  cat("\n")
  
  
  
  vector_of_unique_positions<-unique(input_chunk$pos_for_LDlinkR)
  
  cat("vector_of_unique_positions\n")
  cat(str(vector_of_unique_positions))
  cat("\n")
  
  
  #### subset for what has been already run ----
  
  
  indx.dep<-which(vector_of_unique_positions%in%covered_variants$pos_for_LDlinkR)
  
  if(length(indx.dep) >0){
    
    vector_of_unique_positions_subset<-vector_of_unique_positions[-indx.dep]
    
    
  }else{
    
    vector_of_unique_positions_subset<-vector_of_unique_positions
    
    
  }#length(indx.dep) >0
  
  

  cat("vector_of_unique_positions_subset\n")
  cat(str(vector_of_unique_positions_subset))
  cat("\n")

  #### Run LDlinkR----
  
  
  setwd(path_LDlinkR)
  
  
  result<-LDproxy_batch(snp=vector_of_unique_positions_subset, 
                        pop = "EUR",
                        r2d = "r2",
                        token="1e04f14255a6",
                        genome_build = "grch38",
                        win_size = "500000")
 
}




printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}


#### main script ----

main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--input_chunk"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--covered_variants"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --TranscriptEXP FILE.txt
                        --cadd FILE.txt
                        --ncboost FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  data_wrangling(opt)
  
  
}


###########################################################################

system.time( main() )