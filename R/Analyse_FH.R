#' Analysing FoundHaplo results
#'
#' @description
#' This function analyse the FoundHaplo scores and predict samples that carry the known disease haplotypes for each of the disease-causing variants.
#' @param path_results  Path to a single .txt file with FH scores (type \code{"character"})
#' @param path_to_save_FH_output Path to save the analysis of the FH scores (type \code{"numeric"})
#' @param critical_percentile Critical percentile of the control cohort to derive predictions. Recommend above 99.9 for large cohorts like UKBB. (type \code{"numeric"})
#' @return A dataframe with predicted samples and graphical output will be saved in path_to_save_FH_output

#' \enumerate{
#' \item DCV, name of the disease-causing variant (type \code{"character"})
#' \item test_name, name of the test cohort (type \code{"character"})
#' \item test_sample_name, test sample ID (type \code{"character"})
#' \item FH, corresponding FH score(type \code{"numeric"})
#' }
#' @import dplyr
#' @import ggplot2
#' @import gridExtra
#' @export
#' @examples
#' orig_dir <- getwd()
#' setwd(tempdir())
#' write.table(FH_IBD_scores,paste0(tempdir(),"/results",".txt"),sep = "\t",quote=FALSE, row.names=FALSE,col.names = FALSE) # save FH_IBD_scores
#' Analyse_FH(path_results=paste0(tempdir(),"/results",".txt"),path_to_save_FH_output=tempdir(),critical_percentile=0.99)
#' setwd(orig_dir)

Analyse_FH=function(path_results,path_to_save_FH_output,critical_percentile)
{
  
  results_file=read.delim(path_results,header = FALSE)
  
  colnames(results_file)=c("data_type","test_name","frequency_type","minor_allele_cutoff","imputation_quality_score_cutoff_test",'DCV',"disease_haplotype","test_control_individual","FH_score","left_LLR","right_LLR","total_cM_sharing","total_left_cM_sharing","total_right_cM_sharing","number_of_allele_mismatches_in_the_markov_chain","number_of_markers_in_the_markov_chain","numer_of_haplotype_switches_in_the_markov_chain","snp_density_in_data_file","total_number_of_markers_in_data_file","total_cM_span_of_data_file")
  
  x1=interaction(results_file$DCV,results_file$data_type,results_file$test_name,results_file$disease_haplotype,results_file$test_control_individual)
  results_file$interaction=x1
  
  results_file=distinct(results_file,interaction, .keep_all= TRUE)
  
  
  p1=list()
  p2=list()
  predictions= data.frame(matrix(ncol = 5, nrow = 0))
  
  
  for(ii in 1:length(unique(results_file$DCV)))
  {
    file1=subset(results_file,results_file$DCV %in%  unique(results_file$DCV)[ii])
    file1=file1 %>%
      group_by(DCV,test_name,test_control_individual,data_type) %>%
      summarize(FH = max(FH_score))
    
    file1=as.data.frame(file1)
    file1=na.omit(file1)
    
    attach(file1)
    file1 <- file1[order(-file1$FH),]
    detach(file1)
    
    controls=subset(file1,file1$data_type=="controls")
    test_set=subset(file1,file1$data_type=="test")
    
    
    CLLR=controls %>%
      group_by(test_name) %>%
      summarize(critical_quantile = quantile(FH, probs = critical_percentile)) # you can change the critical value
    
    CLLR=as.data.frame(CLLR)
    
    CLLR=data.frame( x= 1:nrow(CLLR), y = CLLR$critical_quantile)
    
    set.seed(1)
    p1[[ii]]<- ggplot(file1)+
      aes(x=test_name,y=FH,color=data_type)+
      geom_violin(alpha=0.3,position="identity",size=1.5)+scale_color_manual(values=c("black", "salmon3"),labels = c("1000G controls", "Test cohort"))+geom_segment(data = CLLR, color = "black", aes(x = as.numeric(x) - 0.25,y = y,xend = as.numeric(x) + 0.25,yend = y))
    
    p1[[ii]]=  p1[[ii]]+geom_jitter(height=0,width=0.2,alpha=0.5,data = subset(test_set,test_set$FH > CLLR$y),
                                    aes(x=test_name, y=FH),
                                    color = 'red',cex=4,shape=20)
    
    p1[[ii]]=  p1[[ii]]+ ggtitle(paste0("FH score for ",unique(results_file$DCV)[ii]))+theme(plot.title = element_text(color="black", size=14, face="bold"),legend.key.size = unit(0.1, 'cm'),legend.text=element_text(size=10),legend.position="bottom",legend.title=element_blank(),panel.background = element_rect(fill = "white", colour = "white"),axis.text.y=element_text(size=10))+scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+ylab("FH score")
    
    p1[[ii]]=  p1[[ii]]+annotate(geom="text", x=1.4, y=CLLR$y+0.1, label=paste0(critical_percentile," critical percentile"),color="black",size=3)
    
    test_set=subset(file1,file1$data_type=="test")
    test_set=test_set[,-4]
    predicted_samples=subset(test_set,test_set$FH>CLLR$y)
    predictions=data.frame(rbind(predictions,predicted_samples))
    colnames(predictions)=c("DCV","test_name","test_sample_name","FH")
  
    test_set$sample_index=1:nrow(test_set)
    p2[[ii]]=ggplot(test_set)+
      aes(x=sample_index,y=FH)+
      geom_point(color="black")+labs(x="sample index",y="FH score",title="FH score of the test cohort in descending order",subtitle="Alternatively samples sharing the disease haplotypes can be identified if they are seprated as a cluster from the rest")+theme(plot.title = element_text(color="black", size=14, face="bold"),legend.key.size = unit(0.1, 'cm'),legend.text=element_text(size=10),legend.position="bottom",legend.title=element_blank(),panel.background = element_rect(fill = "white", colour = "white"),axis.text.y=element_text(size=10))
    
  }
  
  pdf(paste0(path_to_save_FH_output,"/FH_results.pdf"), height=8, width=14)
  for(ii in 1:length(lengths(p1)))
  {
    print(grid.arrange(p1[[ii]],nrow = 1,ncol=1))
    print(grid.arrange(p2[[ii]],nrow = 1,ncol=1))
    
  }
  dev.off()
  print(paste0("Graphical output is saved in ",path_to_save_FH_output))
  return(predictions)

}

