
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
suppressMessages(library("splitstackshape", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))

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
  
  #### READ and transform LD_thresholds ----
  
  LD_thresholds = as.numeric(unlist(strsplit(opt$LD_thresholds, split=",")))
  
  cat("LD_thresholds_\n")
  cat(sprintf(as.character(LD_thresholds)))
  cat("\n")
  
  #### READ input_file ----
  
  input_file<-as.data.frame(fread(file=opt$input_file, sep="\t", header=T), stringsAsFactors=F)
  
  
  # colnames(input_file)[which(colnames(input_file) == 'Gene Symbol')]<-'Symbol'
  # 
  # input_file$chr<-paste('chr',input_file$chr,sep='')
  
  cat("input_file_0\n")
  cat(str(input_file))
  cat("\n")
  cat(str(unique(input_file$VAR)))
  cat("\n")
  
  
  input_file$chr<-factor(input_file$chr,
                         levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
                                  "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21",
                                  "chr22","chr23","chrX","chrY"), ordered=T)
  
  
  input_file<-input_file[order(input_file$chr,input_file$pos),]

  
  cat("input_file_2\n")
  cat(str(input_file))
  cat("\n")
  cat(str(unique(input_file$chr)))
  cat("\n")
  cat(str(unique(input_file$pos)))
  cat("\n")
  cat(str(unique(input_file$ref)))
  cat("\n")
  cat(str(unique(input_file$alt)))
  cat("\n")
  
  #### Read LDlinkR results ----
  
  path_LDlinkR = paste(out,'LDlinkR','/',sep='')
  
  if (file.exists(path_LDlinkR)){
    
    
  }else{
    
    dir.create(file.path(path_LDlinkR))
    
  }#path_LDlinkR
  
  
  setwd(path_LDlinkR)
  
  files_df<-list.files()
  
  cat("files_df_0\n")
  cat(str(files_df))
  cat("\n")
  
  files_df_sel<-files_df[grep("_grch38\\.txt$",files_df)]
  
  cat("files_df_sel_0\n")
  cat(str(files_df_sel))
  cat("\n")
  
  
  
  
  df_files <- data.frame(matrix(NA,
                                nrow = length(files_df_sel),
                                ncol = 2))
  colnames(df_files)[which(colnames(df_files) == 'X1')]<-'pos_for_LDlinkR'
  colnames(df_files)[which(colnames(df_files) == 'X2')]<-'path'
  
  cat("df_files_0\n")
  cat(str(df_files))
  cat("\n")
  
  
  df_files$path<-paste(path_LDlinkR,files_df_sel,sep='')
  df_files$pos_for_LDlinkR<-gsub("_grch38\\.txt$","",files_df_sel)
  
  
  cat("df_files_1\n")
  cat(str(df_files))
  cat("\n")
  
  
  ##### Reading loop -----
  
  pos_for_LDlinkR_array<-df_files$pos_for_LDlinkR
  
  
  list_results<-list()
  
  list_results_autologous<-list()
  
  
  DEBUG<-0
  
  for(i in 1:length(pos_for_LDlinkR_array))
  {
    pos_for_LDlinkR_array_sel<-pos_for_LDlinkR_array[i]
    
    cat("------------------------------------------------>\t")
    cat(sprintf(as.character(i)))
    cat("\t")
    cat(sprintf(as.character(pos_for_LDlinkR_array_sel)))
    cat("\n")
    
    df_files_sel<-df_files[which(df_files$pos_for_LDlinkR == pos_for_LDlinkR_array_sel),]
    
    if(DEBUG == 1)
    {
      cat("df_files_sel_0\n")
      cat(str(df_files_sel))
      cat("\n")
      
    }
    
    
    
    if(dim(df_files_sel)[1] >0)
    {
      sel_file<-df_files_sel$path
      
      # cat("--->\t")
      # cat(sprintf(as.character(sel_file)))
      # cat("\n")
      
      SIZE_gate<-file.info(sel_file)$size
      
      if(DEBUG == 1)
      {
        cat("SIZE_gate_0\n")
        cat(str(SIZE_gate))
        cat("\n")
      }
      
      if(SIZE_gate> 0)
      {
        LINE_gate<-length(readLines(sel_file))
        
        if(DEBUG == 1)
        {
          cat("LINE_gate_0\n")
          cat(str(LINE_gate))
          cat("\n")
        }
        
        if(LINE_gate> 0)
        {
          results<-as.data.frame(fread(file=sel_file, sep="\t", header = F, skip=1, fill=TRUE), stringsAsFactors = F)
          
          colnames(results)<-c("V1","RS_Number","Coord","Alleles","MAF","Distance","Dprime","R2","Correlated_Alleles","FORGEdb","RegulomeDB","Function")
          
          results$pos_for_LDlinkR<-pos_for_LDlinkR_array_sel
          
          if(DEBUG == 1)
          {
            cat("results_0\n")
            cat(str(results))
            cat("\n")
          }
          
          results<-merge(results,
                         input_file,
                         by='pos_for_LDlinkR')
          
          if(DEBUG == 1)
          {
            cat("results_1\n")
            cat(str(results))
            cat("\n")
          }
          
          
          results_subset<-results[-which(results$RS_Number == results$rs),]
          
          
          if(DEBUG == 1)
          {
            cat("results_subset_1\n")
            cat(str(results_subset))
            cat("\n")
          }
          
          results_autologous<-results[which(results$RS_Number == results$rs),]
          
          if(DEBUG == 1)
          {
            cat("results_autologous_1\n")
            cat(str(results_autologous))
            cat("\n")
          }
          
          list_results_autologous[[i]]<-results_autologous
          
          
          list_thresholds<-list()
          
          for(k in 1:length(LD_thresholds)){
            
            LD_thresholds_sel<-LD_thresholds[k]
            
            # cat("------>\t")
            # cat(sprintf(as.character(LD_thresholds_sel)))
            # cat("\n")
            
            
            results_subset_thresholded<-results_subset[which(results_subset$R2 >= LD_thresholds_sel),]
            
            if(dim(results_subset_thresholded)[1] >0){
              
              results_subset_thresholded$LD_threshold<-LD_thresholds_sel
              
              if(DEBUG == 1)
              {
                cat("results_subset_thresholded_0\n")
                cat(str(results_subset_thresholded))
                cat("\n")
              }
              
              indx.int<-c(which(colnames(results_subset_thresholded) == 'VAR'),which(colnames(results_subset_thresholded) == 'rs'),
                          which(colnames(results_subset_thresholded) == 'Symbol'),
                          which(colnames(results_subset_thresholded) == 'R2'),which(colnames(results_subset_thresholded) == 'RS_Number'),which(colnames(results_subset_thresholded) == 'MAF'),
                          which(colnames(results_subset_thresholded) == 'LD_threshold'))
              
              
              results_subset_thresholded_subset<-unique(results_subset_thresholded[,indx.int])
              
              if(DEBUG == 1)
              {
                cat("results_subset_thresholded_subset_0\n")
                cat(str(results_subset_thresholded_subset))
                cat("\n")
              }
              
              list_thresholds[[k]]<-results_subset_thresholded_subset
              
              
            }#dim(results_subset_thresholded)[1] >0
            
          }# k in 1:length(LD_thresholds)
         
          if(length(list_thresholds) >0)
          {
            results_thresholds = unique(as.data.frame(data.table::rbindlist(list_thresholds, fill = T)))
            
            if(DEBUG == 1)
            {
              cat("results_thresholds_0\n")
              cat(str(results_thresholds))
              cat("\n")
            }
            
            results_thresholds.dt<-data.table(results_thresholds, key=c("VAR","rs","Symbol","LD_threshold"))
            
            
            results_thresholds_collapsed<-as.data.frame(results_thresholds.dt[,.(Proxy_rsid=paste(RS_Number, collapse=";"),
                                                                                  Proxy_R2= paste(R2, collapse=";"),
                                                                                 Proxy_MAF=paste(MAF, collapse=";")), by=key(results_thresholds.dt)],
                                                        stringAsFactors=F)
            
            if(DEBUG == 1)
            {
              cat("results_thresholds_collapsed_0\n")
              cat(str(results_thresholds_collapsed))
              cat("\n")
            }
            
            
            list_results[[i]]<-results_thresholds_collapsed
            
          }#list_thresholds
          
          
        }#LINE_gate> 0
      }# SIZE_gate> 0
    }# dim(df_files_sel)[1] >0
    
  }# i in 1:length(pos_for_LDlinkR_array)
  
  
  if(length(list_results) >0)
  {
    Results = unique(as.data.frame(data.table::rbindlist(list_results, fill = T)))
    
    cat("Results_0\n")
    cat(str(Results))
    cat("\n")
    cat(str(unique(Results$VAR)))
    cat("\n")
    cat(str(unique(Results$rs)))
    cat("\n")
    
    
    # ####chek AF-----
    # 
    # Results_long<-unique(as.data.frame(cSplit(Results,sep = ';', direction = "long",
    #                                                                            splitCols = c("Proxy_rsid","Proxy_R2","Proxy_MAF")),stringsAsFactors=F))
    # 
    # 
    # cat("Results_long_0\n")
    # cat(str(Results_long))
    # cat("\n")
    # 
    # Results_long.dt<-data.table(Results_long, key=c('LD_threshold'))
    # 
    # 
    # 
    # Results_long_MIN<-as.data.frame(Results_long.dt[,.SD[which.min(Proxy_MAF)], by=key(Results_long.dt)], stringsAsFactors=F)
    # 
    # cat("Results_long_MIN_0\n")
    # cat(str(Results_long_MIN))
    # cat("\n")
    
    #### Rescue variants not in 1000G pannel ----
    
    
    input_file_rescue<-input_file[-which(input_file$VAR%in%Results$VAR),]
    
    cat("input_file_rescue_0\n")
    cat(str(input_file_rescue))
    cat("\n")
    cat(str(unique(input_file_rescue$VAR)))
    cat("\n")
    cat(str(unique(input_file_rescue$rs)))
    cat("\n")
    
    
    Results<-merge(Results,
                   input_file_rescue,
                   by=colnames(input_file_rescue)[which(colnames(input_file_rescue)%in%colnames(Results))],
                   all=T)
    
    
    cat("Results_1\n")
    cat(str(Results))
    cat("\n")
    cat(str(unique(Results$VAR)))
    cat("\n")
    cat(str(unique(Results$rs)))
    cat("\n")
    
    
    
   ### save ----
    
    setwd(out)
    
    write.table(Results, file='sQTLs_Proxies.tsv', sep="\t",quote=F,row.names=F)
    
    # setwd(out)
    # 
    # write.table(Results_long_MIN, file='sQTLs_Proxies_Min_MAF.tsv', sep="\t",quote=F,row.names=F)
    
    
  }# length(list_results) >0
  
  
  if(length(list_results_autologous) >0)
  {
    Results = unique(as.data.frame(data.table::rbindlist(list_results_autologous, fill = T)))
    
    cat("Results_0\n")
    cat(str(Results))
    cat("\n")
    cat(str(unique(Results$VAR)))
    cat("\n")
    cat(str(unique(Results$rs)))
    cat("\n")
    
    
    Results.dt<-data.table(Results, key=c('Symbol'))



    Results_MIN<-as.data.frame(Results.dt[,.SD[which.min(MAF)], by=key(Results.dt)], stringsAsFactors=F)

    cat("Results_MIN_0\n")
    cat(str(Results_MIN))
    cat("\n")
    
    ### save ----
    
    setwd(out)
    
    write.table(Results_MIN, file='sQTLs_min_MAF.tsv', sep="\t",quote=F,row.names=F)
    
  }# length(list_results_autologous) >0
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
    make_option(c("--LD_thresholds"), type="character", default=NULL, 
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