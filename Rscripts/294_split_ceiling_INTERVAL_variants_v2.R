suppressMessages(library("plyr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("data.table", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("crayon", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggplot2", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("farver", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("labeling", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("optparse", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("dplyr", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("backports", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("broom", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rstudioapi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cli", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tzdb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggeasy", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("sandwich", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("digest", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tidyverse", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("RColorBrewer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cowplot", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))




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
  
  
  #### READ input_file ----
  
  input_file<-as.data.frame(fread(file=opt$input_file, sep="\t", header=T), stringsAsFactors=F)
  

  cat("input_file_0\n")
  cat(str(input_file))
  cat("\n")
  cat(str(unique(input_file$VAR)))
  cat("\n")
  cat(str(unique(input_file$pos_for_LDlinkR)))
  cat("\n")
  
  #### READ covered_variants ----
  
  covered_variants<-as.data.frame(fread(file=opt$covered_variants, sep="\t", header=F), stringsAsFactors=F)
  
  colnames(covered_variants)<-'pos_for_LDlinkR'
  
  
  cat("covered_variants_0\n")
  cat(str(covered_variants))
  cat("\n")
  
  
  input_file_subset<-input_file[-which(input_file$pos_for_LDlinkR%in%covered_variants$pos_for_LDlinkR),]
  
  
  cat("input_file_subset_0\n")
  cat(str(input_file_subset))
  cat("\n")
  cat(str(unique(input_file_subset$VAR)))
  cat("\n")
  cat(str(unique(input_file_subset$pos_for_LDlinkR)))
  cat("\n")
  
  
  vector_of_unique_positions<-unique(input_file_subset$pos_for_LDlinkR)
  
  cat("vector_of_unique_positions\n")
  cat(str(vector_of_unique_positions))
  cat("\n")
  

  #### chunk it ----
  

  
  List_of_Vars<-split(vector_of_unique_positions, ceiling(seq_along(vector_of_unique_positions)/1000))
  
  cat("List_of_Vars_\n")
  cat(str(List_of_Vars))
  cat("\n")
 
  
  
  ##### Printing loop ---- 
  
  for(i in 1:length(List_of_Vars))
  {
    name_sel<-paste("chunk",i,sep="_")
    chunk_sel<-List_of_Vars[[i]]
    
    cat("------------------------------->name_sel\n")
    cat(sprintf(as.character(name_sel)))
    cat("\n")
    
    
    filename_1<-paste("List_",name_sel,".txt", sep='')
    
    write.table(chunk_sel,
                file=filename_1, sep="\t", quote=F, row.names = F)
    
  }# i in 1:length(List_of_Vars)
 
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
    make_option(c("--input_file"), type="character", default=NULL, 
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