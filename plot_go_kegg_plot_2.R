#!/usr/bin/Rscript
pars <- commandArgs(trailingOnly = TRUE)
if(length(pars)!=7){
  print("plot_kegg_go_plot_2.R prefix(no.filter* file) pmodule(1:pajust; 0:pvalue) pvalue(0.05) width(10) height(30) plot_kegg_go(1:kegg and go; 0:kegg) topnum(20)")
  quit()
}
#########set parameter###############
#######set directory######
#######inputfile name########
input_name<-pars[1]
###p_module: 1 is p.ajust;0 is pvalue##
p_module<-as.numeric(pars[2])
#######p.ajust/pvalue value#########
p_value <- as.numeric(pars[3])
#######pdf width#########
my_width <- as.numeric(pars[4])
#######pdf height#########
my_height <- as.numeric(pars[5])
#####kegg_go_plot:1 plot kegg and GO picture;0 just  plot kegg picture####
kegg_go_plot<-as.numeric(pars[6])
####top number#######
topNum<-as.numeric(pars[7])

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
  filename<-paste0(input_name,"_",enrich_list[iterm],".no_filter.xlsx")
  df<-read.xlsx(filename,sheet = 1)
  df<-read.xlsx(filename)
  df$p.adjust<-as.numeric(df$p.adjust)
  df$pvalue<-as.numeric(df$pvalue)
  getRatio<-function(x){ 
    num<-unlist(strsplit(x,split = "|",fixed = T))
    ratio<-as.numeric(num[1])/as.numeric(num[2])
    ratio<-round(ratio,2)
    return(ratio)
  }
  
  if(p_module==1){
    my_df <- subset(df, select = c(Cluster, ID,Description, p.adjust, Count,GeneRatio))
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
  }else if(p_module==0){
    print(paste0(input_name,"_",enrich_list[iterm],".no_filter",".xlsx"," filtered by pvalue have ",dim(outputfile)[1]))
  }
  
  if(dim(outputfile)[1]==0){
  
    next()
  }else if(dim(outputfile)[1]>0){
    
    write.xlsx(outputfile,paste0(paste0(input_name,"_",enrich_list[iterm]),".filter",p_value,".xlsx"))
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
    }else if(p_module==0){
      print(paste0(input_name,"_",enrich_list[iterm],".no_filter",".xlsx"," filtered by pvalue have not diff pathway/go term",enrich_list[iterm]))
    }
    next()
  }else{
    out_plot_name <- paste0(paste0(input_name,"_",enrich_list[iterm]), '.filter', p_value, '.all.pdf')
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
    out_name <- paste0(paste0(input_name,"_",enrich_list[iterm]), '.filter', p_value, '.top',topNum,'.pdf')
    
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


