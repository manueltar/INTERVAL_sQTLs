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
  
  input_file$variant_chr<-paste('chr',input_file$variant_chr,sep='')
  
  cat("input_file_0\n")
  cat(str(input_file))
  cat("\n")
  cat(str(unique(input_file$VAR)))
  cat("\n")
  
  
  input_file$variant_chr<-factor(input_file$variant_chr,
                       levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
                                "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21",
                                "chr22","chr23","chrX","chrY"), ordered=T)
  
  
  input_file<-input_file[order(input_file$variant_chr,input_file$variant_pos),]
  
  cat("input_file_2\n")
  cat(str(input_file))
  cat("\n")
  cat(str(unique(input_file$variant_chr)))
  cat("\n")
  cat(str(unique(input_file$variant_pos)))
  cat("\n")
  cat(str(unique(input_file$ref)))
  cat("\n")
  cat(str(unique(input_file$alt)))
  cat("\n")
  
  
  colnames(input_file)[which(colnames(input_file) == 'variant_rsid')]<-'rs'
  colnames(input_file)[which(colnames(input_file) == 'variant_chr')]<-'chr'
  colnames(input_file)[which(colnames(input_file) == 'variant_pos')]<-'pos'
  colnames(input_file)[which(colnames(input_file) == 'variant_ref')]<-'ref'
  colnames(input_file)[which(colnames(input_file) == 'variant_alt')]<-'alt'
  colnames(input_file)[which(colnames(input_file) == 'gene_name')]<-'Symbol'
  
  
  
  cat("input_file_3\n")
  cat(str(input_file))
  cat("\n")
  
  
  
  
  indx.int<-c(which(colnames(input_file) == 'chr'),which(colnames(input_file) == 'pos'),
              which(colnames(input_file) == 'ref'),which(colnames(input_file) == 'alt'),
              which(colnames(input_file) == 'Symbol'),which(colnames(input_file) == 'rs'))
  
  input_file_subset<-unique(input_file[,indx.int])
  
  cat("input_file_subset_0\n")
  cat(str(input_file_subset))
  cat("\n")
  
  
  input_file_subset$pos_for_LDlinkR<-paste(input_file_subset$chr,input_file_subset$pos,sep=':')
  input_file_subset$Annotation<-'INTERVAL_sQTL'
  input_file_subset$VAR<-paste(input_file_subset$chr,input_file_subset$pos,input_file_subset$ref,input_file_subset$alt, sep="_")
  
  cat("input_file_subset_1\n")
  cat(str(input_file_subset))
  cat("\n")
  
 
  
  ###### SAVE -----
  
  setwd(out)
  
  write.table(input_file_subset, file="INTERVAL_sQTL.tsv",sep="\t", quote=F, row.names = F)
  saveRDS(input_file_subset, file="INTERVAL_sQTL.rds")
  
  
 
 
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
    make_option(c("--input_file"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--VAR_array"), type="character", default=NULL, 
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