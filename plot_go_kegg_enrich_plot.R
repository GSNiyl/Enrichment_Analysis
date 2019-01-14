#!/usr/bin/Rscript
pars <- commandArgs(trailingOnly = TRUE)
if(length(pars)!=6){
  print("plot_kegg_go_enrich_plot.R inputName outdir pmodule(1:pajust; 0:pvalue) pvalue(0.05) width(10) height(30)")
  quit()
}

#########set parameter###############
########input name#####
file_name<-pars[1]
###output directory##########
out_dir<-pars[2]

###p_module: 1 is p.ajust;0 is pvalue##
p_module<-as.numeric(pars[3])
#######p.ajust/pvalue value#########
p_value <- as.numeric(pars[4])
#######pdf width#########
my_width <- as.numeric(pars[5])
#######pdf height#########
my_height <- as.numeric(pars[6])
#####kegg_go_plot:1 plot kegg and GO picture;0 just  plot kegg picture####
kegg_go_plot<-1
####top number#######
topNum<-20

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(openxlsx)

#####input excel file########
file<-unique(read.table(file_name,sep="\t",header=T,stringsAsFactors=F))
file<-file[is.na(file[,2])==F,]

file_list<-split(file,file[,1])
transf_fun<-function(each_list){
  gene_list1 <- bitr(each_list[,2], fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db', drop = T)$ENTREZID
  return(gene_list1)
}
samples <- lapply(file_list, transf_fun)

###########solve problem which output excel file is bigger ############

##########################KEGG Enrichment analysis#####################

kegg <- compareCluster(samples, fun = 'enrichKEGG', organism = 'hsa', pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = 'BH')
sampleName<-names(samples)
if(length(sampleName)==1){
  foutkegg_xlsx<-paste0(sampleName,"_kegg.no_filter.xlsx")
}else{
  foutkegg_xlsx <- paste0(paste(sampleName,collapse= "_vs_"), '_kegg.no_filter.xlsx')
}


for(i in 1:(nrow(attributes(kegg)$compareClusterResult))){
  attributes(kegg)$compareClusterResult$geneID[i] <- paste(bitr(strsplit(summary(kegg)$geneID[i], '/')[[1]], fromType = 'ENTREZID', toType = 'SYMBOL', OrgDb = 'org.Hs.eg.db', drop = T)$SYMBOL, collapse = '|')
}

attributes(kegg)$compareClusterResult$GeneRatio <- gsub('/', '|', attributes(kegg)$compareClusterResult$GeneRatio)
attributes(kegg)$compareClusterResult$BgRatio <- gsub('/', '|', attributes(kegg)$compareClusterResult$BgRatio)
kegg<-as.data.frame(kegg)
write.xlsx(kegg, file = paste0(out_dir,"/",foutkegg_xlsx) )

##########################GO CC Enrichment analysis##################
go_cc <- compareCluster(samples, fun = 'enrichGO', ont = 'CC', OrgDb = 'org.Hs.eg.db', pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = 'BH')
sampleName<-names(samples)
if(length(sampleName)==1){
  foutgo_cc_xlsx<-paste0(sampleName,"_go_cc.no_filter.xlsx")
}else{
  foutgo_cc_xlsx <- paste0(paste(sampleName,collapse= "_vs_"), '_go_cc.no_filter.xlsx')
}

for(i in 1:(nrow(attributes(go_cc)$compareClusterResult))){
  attributes(go_cc)$compareClusterResult$geneID[i] <- paste(bitr(strsplit(summary(go_cc)$geneID[i], '/')[[1]], fromType = 'ENTREZID', toType = 'SYMBOL', OrgDb = 'org.Hs.eg.db', drop = T)$SYMBOL, collapse = '|')
}

attributes(go_cc)$compareClusterResult$GeneRatio <- gsub('/', '|', attributes(go_cc)$compareClusterResult$GeneRatio)
attributes(go_cc)$compareClusterResult$BgRatio <- gsub('/', '|', attributes(go_cc)$compareClusterResult$BgRatio)
go_cc<-as.data.frame(go_cc)
write.xlsx(as.data.frame(go_cc)[,1:10], file = paste0(out_dir,"/",foutgo_cc_xlsx))


##########################GO BP Enrichment analysis##################
go_bp <- compareCluster(samples, fun = 'enrichGO', ont = 'BP', OrgDb = 'org.Hs.eg.db', pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = 'BH')

sampleName<-names(samples)
if(length(sampleName)==1){
  foutgo_bp_xlsx<-paste0(sampleName,"_go_bp.no_filter.xlsx")
}else{
  foutgo_bp_xlsx <- paste0(paste(sampleName,collapse= "_vs_"), '_go_bp.no_filter.xlsx')
}

for(i in 1:(nrow(attributes(go_bp)$compareClusterResult))){
  attributes(go_bp)$compareClusterResult$geneID[i] <- paste(bitr(strsplit(summary(go_bp)$geneID[i], '/')[[1]], fromType = 'ENTREZID', toType = 'SYMBOL', OrgDb = 'org.Hs.eg.db', drop = T)$SYMBOL, collapse = '|')
}

attributes(go_bp)$compareClusterResult$GeneRatio <- gsub('/', '|', attributes(go_bp)$compareClusterResult$GeneRatio)
attributes(go_bp)$compareClusterResult$BgRatio <- gsub('/', '|', attributes(go_bp)$compareClusterResult$BgRatio)
go_bp<-as.data.frame(go_bp)
write.xlsx(as.data.frame(go_bp)[,1:10], file = paste0(out_dir,"/",foutgo_bp_xlsx))

##########################GO MF Enrichment analysis##################
go_mf <- compareCluster(samples, fun = 'enrichGO', ont = 'MF', OrgDb = 'org.Hs.eg.db', pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = 'BH')

sampleName<-names(samples)
if(length(sampleName)==1){
  foutgo_mf_xlsx<-paste0(sampleName,"_go_mf.no_filter.xlsx")
}else{
  foutgo_mf_xlsx <- paste0(paste(sampleName,collapse= "_vs_"), '_go_mf.no_filter.xlsx')
}
for(i in 1:(nrow(attributes(go_mf)$compareClusterResult))){
  attributes(go_mf)$compareClusterResult$geneID[i] <- paste(bitr(strsplit(summary(go_mf)$geneID[i], '/')[[1]], fromType = 'ENTREZID', toType = 'SYMBOL', OrgDb = 'org.Hs.eg.db', drop = T)$SYMBOL, collapse = '|')
}

attributes(go_mf)$compareClusterResult$GeneRatio <- gsub('/', '|', attributes(go_mf)$compareClusterResult$GeneRatio)
attributes(go_mf)$compareClusterResult$BgRatio <- gsub('/', '|', attributes(go_mf)$compareClusterResult$BgRatio)
go_mf<-as.data.frame(go_mf)
write.xlsx(as.data.frame(go_mf)[,1:10], file = paste0(out_dir,"/",foutgo_mf_xlsx))

#######################output enrichment Results#################################


print(paste0("Numbers 0f KEGG enrichment:",dim(kegg)[1]))
print(paste0("Numbers 0f GO_BP enrichment:",dim(go_bp)[1]))
print(paste0("Numbers 0f GO_MF enrichment:",dim(go_mf)[1]))
print(paste0("Numbers 0f GO_MF enrichment:",dim(go_mf)[1]))


a<-paste0("Numbers 0f KEGG enrichment:",dim(kegg)[1])
b<-paste0("Numbers 0f GO_BP enrichment:",dim(go_bp)[1])
c<-paste0("Numbers 0f GO_MF enrichment:",dim(go_mf)[1])
d<-paste0("Numbers 0f GO_MF enrichment:",dim(go_mf)[1])

logfile<-rbind(a,b,c,d)

input_name<-paste(sampleName,collapse= "_vs_")


##############################plot##################################################
########library packages#######
library(openxlsx)
library(stringr) 
library(ggplot2)
if(kegg_go_plot==1){
  enrich_list<-c("kegg","go_bp","go_cc","go_mf")
}else if(kegg_go_plot==0){
  enrich_list<-c("kegg")
}


for(iterm in c(1:length(enrich_list))){
  if(enrich_list[iterm]=="kegg"){
    df<-kegg
  }else if(enrich_list[iterm]=="go_bp"){
    df<-go_bp
  }else if(enrich_list[iterm]=="go_cc"){
    df<-go_cc
  }else if(enrich_list[iterm]=="go_mf"){
    df<-go_mf
  }
  
  getRatio<-function(x){ 
    num<-unlist(strsplit(x,split = "|",fixed = T))
    ratio<-as.numeric(num[1])/as.numeric(num[2])
    ratio<-round(ratio,2)
    return(ratio)
  }
  
  if(p_module==1){
    my_df <- subset(df, select = c(Cluster,ID,Description,GeneRatio, p.adjust, Count))
    my_df$Cluster<-factor(my_df$Cluster,levels = sort(unique(my_df$Cluster)))
    df.select<-my_df[my_df$p.adjust<p_value,]
    
    outputfile<-df[df$p.adjust<=p_value,]
  }else if(p_module==0){
    my_df <-  subset(df, select = c(Cluster, ID,Description, pvalue, Count,GeneRatio))
    my_df$Cluster<-factor(my_df$Cluster,levels = sort(unique(my_df$Cluster)))
    df.select<-my_df[my_df$pvalue<p_value,]
    outputfile<-df[df$pvalue<=p_value,]
  }
  
  if(p_module==1){
    print(paste0(input_name,"_",enrich_list[iterm],".no_filter",".xlsx"," filtered by p.ajust have ",dim(outputfile)[1]))
    
    aa<-paste0(input_name,"_",enrich_list[iterm],".no_filter",".xlsx"," filtered by p.ajust have ",dim(outputfile)[1])
  }else if(p_module==0){
    print(paste0(input_name,"_",enrich_list[iterm],".no_filter",".xlsx"," filtered by pvalue have ",dim(outputfile)[1]))
    aa<-paste0(input_name,"_",enrich_list[iterm],".no_filter",".xlsx"," filtered by pvalue have ",dim(outputfile)[1])
  }
  logfile<-rbind(logfile,aa)
  
  if(dim(outputfile)[1]==0){

    next()
  }else if(dim(outputfile)[1]>0){
    
    write.xlsx(outputfile,paste0(paste0(out_dir,"/",input_name,"_",enrich_list[iterm]),".filter",p_value,".xlsx"))
  }
  
  df.select$GeneRatio<-as.numeric(as.matrix(lapply(df.select$GeneRatio,getRatio)))
  df.select<-df.select[order(df.select$Description),]  
  sample.num<-unique(df$Cluster)
  
  if(length(sample.num)<=2){
    deduplicate.pathway.df<-df.select[!(duplicated(df.select$ID)|duplicated(df.select$ID, fromLast = TRUE)),]
  }else{
    uniq.id<-unique(df.select$ID)
    for(i in 1:length(uniq.id)){
      delete.duplicate<-df.select[df.select$ID%in%df.select$ID[i],]
      if(dim(delete.duplicate)[1]==length(sample.num)){
        df.select<-df.select[!df.select$ID%in%df.select$ID[i],]
      }
      
    }
    deduplicate.pathway.df<-df.select
    
  }
  deduplicate.pathway.df$Count <- as.integer(deduplicate.pathway.df$Count)
  if(p_module==1){
    deduplicate.pathway.df<-deduplicate.pathway.df[order(deduplicate.pathway.df$p.adjust),]
  }else if(p_module==0){
    deduplicate.pathway.df<-deduplicate.pathway.df[order(deduplicate.pathway.df$pvalue),]
  }
  
  n <- nrow(deduplicate.pathway.df)
  if(n <= 10){
    my_height <- 5
  }else if(n <= 50){
    my_height <- 8
  }else if(n>50&&n<=120){
    my_height <- as.numeric(my_height)
  }else if(n>120 & n<1000){
    my_height<-50
  }else if(n > 1000& n<2000){
    my_height<-120
  }else if(n>=2000& n<3000){
    my_height<-200
  }else if(n>3000){
    my_height<-300
    my_width<-13
  }
  
  ########################draw all different pathway##################
  if(dim(deduplicate.pathway.df)[1]==0){
    if(p_module==1){
      print(paste0(input_name,"_",enrich_list[iterm],".no_filter",".xlsx"," filtered by p.ajust have not different pathway/go term",enrich_list[iterm]))
      bb<-paste0(input_name,"_",enrich_list[iterm],".no_filter",".xlsx"," filtered by p.ajust have not diff pathway/go term",enrich_list[iterm])
    }else if(p_module==0){
      print(paste0(input_name,"_",enrich_list[iterm],".no_filter",".xlsx"," filtered by pvalue have not diff pathway/go term",enrich_list[iterm]))
      bb<-paste0(input_name,"_",enrich_list[iterm],".no_filter",".xlsx"," filtered by pvalue have not diff pathway/go term",enrich_list[iterm])
    }
    logfile<-rbind(logfile,aa,bb)
    next()
  }else{
    out_plot_name <- paste0(paste0(out_dir,"/",input_name,"_",enrich_list[iterm]), '.filter', p_value, '.all.pdf')
    if(p_module==1){
      p <- ggplot(deduplicate.pathway.df, aes(x = Cluster, y = Description)) +
        geom_point(aes(size = GeneRatio, color = p.adjust)) +
        scale_colour_gradient(low = 'red', high = 'blue') +
        scale_size(range = c(2, 5)) +
        xlab('') +
        ylab('') +
        labs(size='GeneRatio', colour='p.adjust') +
        scale_x_discrete(drop = F)
      ggsave(filename = out_plot_name, width = as.numeric(my_width), height = as.numeric(my_height), limitsize = F)
    }else if(p_module==0){
      p <- ggplot(deduplicate.pathway.df, aes(x = Cluster, y = Description)) +
        geom_point(aes(size = GeneRatio, color = pvalue)) +
        scale_colour_gradient(low = 'red', high = 'blue') +
        scale_size(range = c(2, 5)) +
        xlab('') +
        ylab('') +
        labs(size='GeneRatio', colour='p.adjust') +
        scale_x_discrete(drop = F)
      ggsave(filename = out_plot_name, width = as.numeric(my_width), height = as.numeric(my_height), limitsize = F)
    }
    
    if(length(sample.num)<=2){
      if(dim(deduplicate.pathway.df)[1]<=topNum){
        deduplicate.pathway.df.top<-deduplicate.pathway.df
      }else{
        deduplicate.pathway.df.top<-deduplicate.pathway.df[1:topNum,]
      }
      
    }else{
      deduplicate.pathway.df1<-deduplicate.pathway.df
      if(dim(deduplicate.pathway.df)[1]<=topNum){
        deduplicate.pathway.df.top<-deduplicate.pathway.df1
      }else{
        df_list <- split(deduplicate.pathway.df1, deduplicate.pathway.df1$Cluster)
        df_list<-lapply(df_list,function(x) x[order(x[,4]),])
        uniq_list<-unlist(lapply(df_list,function(x){dim(x)[1]}))
        subset_sample<-names(uniq_list[uniq_list[]==max(uniq_list)])
        if(length(subset_sample)>1){
          plot_sample_df<-deduplicate.pathway.df1[deduplicate.pathway.df1$Cluster==subset_sample[1],]
        }else{
          plot_sample_df<-deduplicate.pathway.df1[deduplicate.pathway.df1$Cluster==subset_sample,]
        }
        if(dim(plot_sample_df)[1]>topNum){
          plot_sample_df_tmp<-plot_sample_df[1:topNum,]
        }else{
          plot_sample_df_tmp<-plot_sample_df
        }
        my_fun <- function(each_list){
          tmp_df <- each_list[each_list['Count'] != 0,]
          tmp_df<-tmp_df[match(plot_sample_df_tmp$Description,tmp_df$Description,0L),]
          return(tmp_df)
        }
        
        df_list_filter<- lapply(df_list, my_fun)
        deduplicate.pathway.df.top <- do.call('rbind', df_list_filter)
      }
    }
    n <- nrow(deduplicate.pathway.df.top)
    if(n <= 10){
      my_height <- 5
    }else if(n <20){
      my_height <- 5
    }else if(n<30){
      my_height <- 7
    }else if(n <= 50){
      my_height <- 8
    }else{
      my_height<-10
    }
    ########################draw top20 different pathway##################
    out_name <- paste0(paste0(out_dir,"/",input_name,"_",enrich_list[iterm]), '.filter', p_value, '.top',topNum,'.pdf')
    if(p_module==1){
      p <- ggplot(deduplicate.pathway.df.top, aes(x =Cluster, y = Description)) +
        geom_point(aes(size = GeneRatio, color = p.adjust)) +
        scale_colour_gradient(low = 'red', high = 'blue') +
        scale_size(range = c(2, 5)) +
        xlab('') +
        ylab('') +
        labs(size='GeneRatio', colour='p.adjust') +
        scale_x_discrete(drop = F)
      ggsave(filename = out_name, width = as.numeric(my_width), height = as.numeric(my_height), limitsize = F)
    }else if(p_module==0){
      p <- ggplot(deduplicate.pathway.df.top, aes(x =Cluster, y = Description)) +
        geom_point(aes(size = GeneRatio, color = pvalue)) +
        scale_colour_gradient(low = 'red', high = 'blue') +
        scale_size(range = c(2, 5)) +
        xlab('') +
        ylab('') +
        labs(size='GeneRatio', colour='p.adjust') +
        scale_x_discrete(drop = F)
      ggsave(filename = out_name, width = as.numeric(my_width), height = as.numeric(my_height), limitsize = F)
    }
  
  }
}
  write.table(logfile,"workLOG.txt",sep = "\t",col.names = F,row.names = F,quote = F)
