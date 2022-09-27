###Analyzing a single cohort
library(dplyr)
results_path="PATH/results.txt" # Add the path to results file
path_to_save_plots="/wehisan/bioinf/lab_bahlo/users/robertson.e/R_package/Example/FH_results"

results_file=read.delim(results_path,header = FALSE)

colnames(results_file)=c("data_type","test_name","frequency_type","minor_allele_cutoff","imputation_quality_score_cutoff_test",'DCV',"disease_haplotype","test_control_individual","FH_score","left_LLR","right_LLR","total_cM_sharing","total_left_cM_sharing","total_right_cM_sharing","number_of_allele_mismatches_in_the_markov_chain","number_of_markers_in_the_markov_chain","numer_of_haplotype_switches_in_the_markov_chain","snp_density_in_data_file","total_number_of_markers_in_data_file","total_cM_span_of_data_file")

x1=interaction(results_file$DCV,results_file$data_type,results_file$test_name,results_file$disease_haplotype,results_file$test_control_individual)
results_file$interaction=x1

results_file=distinct(results_file,interaction, .keep_all= TRUE)


p1=list()
predictions= data.frame(matrix(ncol = 5, nrow = 0))
critical_percentile=0.99 # user input

for(x in 1:length(unique(results_file$DCV)))
  {
  file1=subset(results_file,results_file$DCV %in%  unique(results_file$DCV)[x])
  file1=file1 %>%
    group_by(DCV,test_name,test_control_individual,data_type) %>%
    summarize(CFH = max(FH_score))

  file1=as.data.frame(file1)
  file1=file1 %>% drop_na()

  controls=subset(file1,file1$data_type=="controls")

  CLLR=controls %>%
    group_by(test_name) %>%
    summarize(critical_quantile = quantile(CFH, probs = critical_percentile)) # you can change the critical value

  CLLR=as.data.frame(CLLR)

  CLLR=data.frame( x= 1:nrow(CLLR), y = CLLR$critical_quantile)

  p1[[x]]<- ggplot(file1)+
    aes(x=test_name,y=CFH,color=data_type)+
    geom_violin(alpha=0.3,position="identity")+scale_color_manual(values=c("black", "salmon3"),labels = c("1000G EUR controls (n=398)", "Test cohort (n=100)"))+geom_segment(data = CLLR, color = "black", aes(x = as.numeric(x) - 0.25,y = y,xend = as.numeric(x) + 0.25,yend = y))

  p1[[x]]=  p1[[x]]+ ggtitle(paste0("Combined FH score for ",unique(results_file$DCV)[x]))+theme(plot.title = element_text(color="black", size=14, face="bold"),legend.key.size = unit(0.1, 'cm'),legend.text=element_text(size=10),legend.position="bottom",legend.title=element_blank(),panel.background = element_rect(fill = "white", colour = "white"),axis.text.y=element_text(size=10))+scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+ylab("CFH score")

  p1[[x]]=  p1[[x]]+annotate(geom="text", x=1.4, y=CLLR$y+0.1, label=paste0(critical_percentile," critical percentile"),color="black",size=3)

  test_set=subset(file1,file1$data_type=="test")
  predicted_samples=subset(test_set,test_set$CFH>CLLR$y)
  predictions=data.frame(rbind(predictions,predicted_samples))



}



pdf(paste0(path_to_save_plots,"/FH_results.pdf"), height=8, width=14)
for(i in 1:length(lengths(p1)))
{
  print(grid.arrange(p1[[i]],nrow = 1,ncol=1))
}
dev.off()
predictions
