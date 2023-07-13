#' Analyse FoundHaplo FH score values for predictions
#'
#' @description
#' Predictions that gave FH scores greater than a selected critical percentile of the 1000 Genomes control cohort will be lislted.
#' ECDF of the FH scores are plotted in save_FH_output_DIR for further analysis
#' @param results_FILE Path to a single .txt file with all the FH scores by FoundHaplo
#' @param save_FH_output_DIR Directory to save the graphical output of the FH scores
#' @param critical_percentile Critical percentile of the control cohort to derive predictions
#' @param from_control If the critical value should be calculated from a control cohort or not. Default is from_control=TRUE. We recommend running a control cohort.
#' @return A data frame with predictions that gave FH scores greater than a selected critical percentile. This function plots the predicted samples in a pdf file in save_FH_output_DIR
#' @import dplyr
#' @import ggplot2
#' @import gridExtra
#' @export
#' @examples
#' orig_DIR <- getwd()
#' temp_DIR <- tempdir()
#' setwd( temp_DIR )
#' write.table(FH_IBD_scores,"FH_IBD_scores.txt",sep = "\t",quote=FALSE, row.names=FALSE,col.names = TRUE)
#' Analyse_FH("FH_IBD_scores.txt",temp_DIR,0.99)
#' setwd(orig_DIR)

Analyse_FH=function(results_FILE,save_FH_output_DIR,critical_percentile=0.99,from_control=TRUE)
{

  results_file=read.delim(results_FILE,header = FALSE)

  colnames(results_file)=c("data_type","test_name","frequency_type","minor_allele_cutoff","imputation_quality_score_cutoff_test",'DCV',"disease_haplotype","test_control_individual","FH_score","left_LLR","right_LLR","total_cM_sharing","total_left_cM_sharing","total_right_cM_sharing","number_of_allele_mismatches_in_the_markov_chain","number_of_markers_in_the_markov_chain","numer_of_haplotype_switches_in_the_markov_chain","snp_density_in_data_file","total_number_of_markers_in_data_file","total_cM_span_of_data_file")

  dummy=interaction(results_file$DCV,results_file$data_type,results_file$test_name,results_file$disease_haplotype,results_file$test_control_individual)
  results_file$interaction=dummy

  results_file=distinct(results_file,interaction, .keep_all= TRUE) # remove duplicate entries

  p1=list() # save violin plots
  p2=list() # save ECDF plots

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

if(from_control){
  print("Critical values are calculated from a control cohort")
  if(nrow(controls)==0){stop("There are no controls to calculate the critical values. Please run a control cohort.")}
  CLLR=controls %>%
    group_by(test_name) %>%
    summarize(critical_quantile = quantile(FH, probs = critical_percentile)) # you can change the critical value


}
if(!from_control){
  print("Critical values are calculated from the test cohort")
  CLLR=test_set %>%
    group_by(test_name) %>%
    summarize(critical_quantile = quantile(FH, probs = critical_percentile)) # you can change the critical value

    }

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

  pdf(paste0(save_FH_output_DIR,"/FH_results_",unique(results_file$DCV)[ii],".pdf"), height=8, width=14)
    print(grid.arrange(p1[[ii]],nrow = 1,ncol=1))

  dev.off()
  }
  print(paste0("Graphical output is saved in ",save_FH_output_DIR))

  print(paste0("Test samples that gave FH score values above the ",critical_percentile*100,"th critical percentile are below"))
  return(predictions)

}
