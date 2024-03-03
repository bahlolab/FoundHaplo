
args = commandArgs(trailingOnly=TRUE)

gen_error=args[1] # 1 if genotype errors should be simulated, 0 otherwise
S_Er=args[2] # 1 if switch errors should be simulated, 0 otherwise
IBD_version=args[3] # 1 for simulating and running FoundHaplo for cases and 1.1 for for simulating and running FoundHaplo for conttrols
gen_error_rate=args[4] # Genotype and imputation error rate allowed; default is 0.01 (1%)
S_Er_rate=args[5] # switch error rate allowed; default is 20.05 switches per mbp
path_to_save=args[6] # Directory path to save the output of FoundHaplo IBD sharing for further analysis
meiosis=args[7] # Estimated number of meiosis between disease-test pair; default is 1

root= args[8]  # The folder to 1000 Genomes VCF files in .RDS format. save all 1000 Genomes VCF files by disease_variant (example: FAME1.chr8.vcf.RDS) in RDS format in a single folder.
path_MIS_ref= args[9] # The folder to MIS-1000 Genomes reference panel, which is used to filter in genetic markers based on 1000 Genomes reference panel in the Michigan Imputation Server. Make sure the reference files are in chr.vcf format for all the chromosomes
path_hapmap= args[10] # Directory path to genetic_map_HapMapII_GRCh37 files, which are also FoundHaplo/input_files/public_data/genetic_map_HapMapII_GRCh37/
path_gnomad_frq= args[11] #The path to .txt file that has gnomAD EUR frequency. path_gnomad_frq would be in "/path/hg19_EUR.sites.2015_08.txt" format (This MAF file was taken from ANNOVAR)


seed_id = args[12] # Integer to identify the founder scenario, values take 1-10
sim_ID = args[13] # Integer to track each simulation (1-1320)
RE_loci = args[14] # Name of the disease-causing variant of interest, i.e. FAME1.chr8.119379052. Use the OMIM abbreviation for the disease. Mostly repeat expansion d
sharing_w = args[15] # Expected sharing length in total around the disease variant (0.5,1,2 and 5 in cM)




class(sharing_w)="numeric"
class(sim_ID)="numeric"
class(gen_error)="numeric"
class(S_Er)="numeric"
class(gen_error_rate)="numeric"
class(S_Er_rate)="numeric"
class(meiosis)="numeric"



source("/path/Simulations_FoundHaplo_single_founder_effects.R")



Simulations_FoundHaplo_single_founder_effects(seed_id,sim_ID,RE_loci,sharing_w,gen_error,S_Er,IBD_version,gen_error_rate,S_Er_rate,path_to_save,meiosis,root,path_MIS_ref,path_hapmap,path_gnomad_frq)




