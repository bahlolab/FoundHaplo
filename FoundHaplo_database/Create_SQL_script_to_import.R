Create_SQL_script_to_import=function(disease_hap_FILE,save_SQL_FILE,db_port,db_host,db_password,db_name,db_unix_socket,family_id,individual_id,father_id,mother_id,sex,sex_method,ancestral_population,ancestral_superpopulation,ancestry_method,sample_id,data_type,external_lab_id,external_source,phasing_method,impute_method,impute_phasing_panel,import_date,DCV_id,disease,disease_name,omim_id,gene,genomic_region,inheritance_model,chromosome,start_position_hg19,end_position_hg19,start_position_hg38,end_position_hg38,start_position_cM,end_position_cM,genotype,validated,validation_method,validation_note)
{
  
  library(data.table)
  library(DBI)
  library(RMariaDB)
  disease_hap_FILE <-fread(disease_hap_FILE,skip = "#CHROM")
  
  options(scipen=99)
  
  individuals <- data.frame(family_id=family_id, individual_id=individual_id, father_id=father_id, mother_id=mother_id,
                            sex=sex,sex_method=paste0("\"",sex_method,"\""),ancestral_population=paste0("\"",ancestral_population,"\""),
                            ancestral_superpopulation=paste0("\"",ancestral_superpopulation,"\""),ancestry_method=paste0("\"",ancestry_method,"\""),
                            stringsAsFactors=FALSE)
  
  
  samples <- data.frame(sample_id=sample_id, individual_id=individual_id,data_type=paste0("\"",data_type,"\""),external_lab_id=paste0("\"",external_lab_id,"\""),external_source=paste0("\"",external_source,"\""),phasing_method=paste0("\"",phasing_method,"\""),impute_method=paste0("\"",impute_method,"\""),impute_phasing_panel=paste0("\"",impute_phasing_panel,"\""),import_date=paste0("\"",import_date,"\""),
                        stringsAsFactors=FALSE)
  
  pathogenic_mutations <- data.frame(DCV_id=DCV_id, disease_name=paste0("\"",disease_name,"\""),
                                     omim_id=omim_id, gene=paste0("\"",gene,"\""), genomic_region=paste0("\"",genomic_region,"\""),inheritance_model=paste0("\"",inheritance_model,"\""),
                                     chromosome=paste0("\"",chromosome,"\""), start_position_hg19=start_position_hg19, end_position_hg19=end_position_hg19,start_position_hg38=start_position_hg38,end_position_hg38=end_position_hg38,
                                     start_position_cM=start_position_cM, end_position_cM=end_position_cM,
                                     stringsAsFactors=FALSE)
  
  individuals_with_known_mutations <- data.frame(individual_id=individual_id, DCV_id=DCV_id, genotype=genotype,
                                                 validated=validated, validation_method=paste0("\"",validation_method,"\""),
                                                 validation_note=paste0("\"",validation_note,"\""),
                                                 stringsAsFactors=FALSE)
  
  R2=sapply(strsplit(disease_hap_FILE$INFO,";",fixed=TRUE),"[[", 1)
  if(sum(R2 %like% "R2")==0){
    print("Disease haplotypes are not imputed, set R-squared to zero")
    R2=rep(0,nrow(disease_hap_FILE))
    
    fields=c("AF_raw", "AF_afr", "AF_sas","AF_amr", "AF_eas", "AF_nfe", "AF_fin") # make sure you only have these fields in the VCF file
    
    get_info_fields <- function(info, fields) {
      info_split <- strsplit(strsplit(info, ";")[[1]], "=")
      info_split_names <- sapply(info_split, "[", 1)
      info_split_names[1:6]=c("AF_raw", "AF_afr", "AF_sas","AF_amr", "AF_eas", "AF_nfe", "AF_fin")
      info_split_values <- sapply(info_split, "[", 2)
      
      return(info_split_values[match(fields, fields)])
    }
    
    hap_info_extract <- do.call(rbind, lapply(disease_hap_FILE$INFO, function(x) get_info_fields(x, fields)))
    hap_info_extract <- apply(hap_info_extract, 2, as.numeric)
    hap_info_extract <- apply(hap_info_extract, 2, function(x) ifelse(is.na(x), 0, x))
    hap_info_extract <- as.data.frame(hap_info_extract, stringsAsFactors=FALSE)
    rownames(hap_info_extract) <- rownames(disease_hap_FILE)
    colnames(hap_info_extract) <- fields
    
    disease_hap_FILE=as.data.frame(disease_hap_FILE)
    
    genetic_markers <- data.frame(marker_id=1:nrow(disease_hap_FILE), rs_id=paste0("\"", disease_hap_FILE$ID, "\""),
                                  chromosome=paste0("\"chr", disease_hap_FILE[,"#CHROM"], "\""),
                                  position_hg19=disease_hap_FILE$POS,
                                  position_hg38=-99999,
                                  position_cM=-99999,
                                  reference_allele=paste0("\"", disease_hap_FILE$REF, "\""),
                                  alternate_allele=paste0("\"", disease_hap_FILE$ALT, "\""),
                                  marker_type="\"SNP\"",
                                  maf_gnomad_ALL=as.numeric(hap_info_extract$AF_raw),
                                  maf_gnomad_AFR=as.numeric(hap_info_extract$AF_afr),
                                  maf_gnomad_NFE=as.numeric(hap_info_extract$AF_nfe),
                                  maf_gnomad_FIN=as.numeric(hap_info_extract$AF_fin),
                                  maf_gnomad_AMR=as.numeric(hap_info_extract$AF_amr),
                                  maf_gnomad_EAS=as.numeric(hap_info_extract$AF_eas),
                                  maf_gnomad_SAS=as.numeric(hap_info_extract$AF_sas),
                                  stringsAsFactors=FALSE)
    
  }
  
  # script is different when R-squared does not exist in the VCF file
  if(sum(R2 %like% "R2")==nrow(disease_hap_FILE)){
    R2=sapply(strsplit(R2,"=",fixed=TRUE),"[[", 2)
    
    fields=c("R2","AF_raw", "AF_afr", "AF_sas","AF_amr", "AF_eas", "AF_nfe", "AF_fin") # make sure you only have these fields in the VCF file
    
    get_info_fields <- function(info, fields) {
      info_split <- strsplit(strsplit(info, ";")[[1]], "=")
      info_split_names <- sapply(info_split, "[", 1)
      info_split_names[1:6]=c("R2","AF_raw", "AF_afr", "AF_sas","AF_amr", "AF_eas", "AF_nfe", "AF_fin") # make sure you only have these fields in the VCF file
      info_split_values <- sapply(info_split, "[", 2)
      
      return(info_split_values[match(fields, fields)])
    }
    
    hap_info_extract <- do.call(rbind, lapply(disease_hap_FILE$INFO, function(x) get_info_fields(x, fields)))
    hap_info_extract <- apply(hap_info_extract, 2, as.numeric)
    hap_info_extract <- apply(hap_info_extract, 2, function(x) ifelse(is.na(x), 0, x))
    hap_info_extract <- as.data.frame(hap_info_extract, stringsAsFactors=FALSE)
    rownames(hap_info_extract) <- rownames(disease_hap_FILE)
    colnames(hap_info_extract) <- fields
    
    disease_hap_FILE=as.data.frame(disease_hap_FILE)
    
    genetic_markers <- data.frame(marker_id=1:nrow(disease_hap_FILE), rs_id=paste0("\"", disease_hap_FILE$ID, "\""),
                                  chromosome=paste0("\"chr", disease_hap_FILE[,"#CHROM"], "\""),
                                  position_hg19=disease_hap_FILE$POS,
                                  position_hg38=-99999,
                                  position_cM=-99999,
                                  reference_allele=paste0("\"", disease_hap_FILE$REF, "\""),
                                  alternate_allele=paste0("\"", disease_hap_FILE$ALT, "\""),
                                  marker_type="\"SNP\"",
                                  maf_gnomad_ALL=as.numeric(hap_info_extract$AF_raw),
                                  maf_gnomad_AFR=as.numeric(hap_info_extract$AF_afr),
                                  maf_gnomad_NFE=as.numeric(hap_info_extract$AF_nfe),
                                  maf_gnomad_FIN=as.numeric(hap_info_extract$AF_fin),
                                  maf_gnomad_AMR=as.numeric(hap_info_extract$AF_amr),
                                  maf_gnomad_EAS=as.numeric(hap_info_extract$AF_eas),
                                  maf_gnomad_SAS=as.numeric(hap_info_extract$AF_sas),
                                  stringsAsFactors=FALSE)
    
  }
  
  genotypes <- data.frame(marker_id=1:nrow(disease_hap_FILE), sample_id=sample_id, genotype=disease_hap_FILE$h1,imputed=1,imputation_quality=R2,
                          stringsAsFactors=FALSE)
  
  genotypes <- genotypes[!is.na(genotypes$genotype), ]
  dummy = strsplit(genotypes$genotype,":")[[1]]
  
  dummy=sapply(strsplit(genotypes$genotype,":",fixed=TRUE),"[[", 1)
  dummy=sapply(strsplit(dummy,"/",fixed=TRUE),"[[", 1)
  
  genotypes$genotype=dummy
  
  mysql_commands <- save_SQL_FILE
  
  # connecting to the database
  
  db = dbConnect(RMariaDB::MariaDB(),bigint = 'integer',port=db_port,host=db_host,user ='remote_usr',password=db_password,dbname=db_name,unix.socket=db_unix_socket)
  
  fetch_genetic_markers=dbSendQuery(db, "SELECT * FROM GeneticMarkers;") # can not add LIMIT here as in SQL
  fetch_genetic_markers <- dbFetch(fetch_genetic_markers,)
  
  ii_command <- paste0("USE ",db_name, ";\n")
  cat(ii_command, file=mysql_commands)
  
  # writing the sql commands to import the first disease haplotype into save_SQL_FILE
  if(nrow(fetch_genetic_markers)==0)
  {
    for (ii in seq_len(nrow(individuals))) {
      ii_command <- paste0("INSERT INTO Individuals (", paste(names(individuals), collapse=","), ") VALUES(", paste(individuals[ii, ], collapse=","), ");\n")
      cat(ii_command, file=mysql_commands, append=TRUE)
    }
    
    for (ii in seq_len(nrow(samples))) {
      ii_command <- paste0("INSERT INTO Samples (", paste(names(samples), collapse=","), ") VALUES(", paste(samples[ii, ], collapse=","), ");\n")
      cat(ii_command, file=mysql_commands, append=TRUE)
    }
    
    for (ii in seq_len(nrow(pathogenic_mutations))) {
      ii_command <- paste0("INSERT INTO DiseaseCausingVariants (", paste(names(pathogenic_mutations), collapse=","), ") VALUES(", paste(pathogenic_mutations[ii, ], collapse=","), ");\n")
      cat(ii_command, file=mysql_commands, append=TRUE)
    }
    
    for (ii in seq_len(nrow(individuals_with_known_mutations))) {
      ii_command <- paste0("INSERT INTO IndividualsWithDiseaseCausingVariants (", paste(names(individuals_with_known_mutations), collapse=","), ") VALUES(", paste(individuals_with_known_mutations[ii, ], collapse=","), ");\n")
      cat(ii_command, file=mysql_commands, append=TRUE)
    }
    
    for (ii in seq_len(nrow(genetic_markers))) {
      ii_command <- paste0("INSERT INTO GeneticMarkers (", paste(names(genetic_markers), collapse=","), ") VALUES(", paste(genetic_markers[ii, ], collapse=","), ");\n")
      cat(ii_command, file=mysql_commands, append=TRUE)
    }
    
    for (ii in seq_len(nrow(genotypes))) {
      ii_command <- paste0("INSERT INTO Genotypes (", paste(names(genotypes), collapse=","), ") VALUES(", paste(genotypes[ii, ], collapse=","), ");\n")
      cat(ii_command, file=mysql_commands, append=TRUE)
    }
  }
  # writing the sql commands to import the first disease haplotype into save_SQL_FILE
  
  # writing the sql commands to import second disease onwards into save_SQL_FILE
  if(nrow(fetch_genetic_markers)>0)
  {
    
    end_marker=max(fetch_genetic_markers$marker_id)
    
    rs_id=paste0("\"", fetch_genetic_markers$rs_id, "\"")
    rs_id=as.data.frame(rs_id,stringsAsFactors=FALSE)
    
    fetch_genetic_markers$rs_id=rs_id
    
    chromosome=paste0("\"", fetch_genetic_markers$chromosome, "\"")
    chromosome=as.data.frame(chromosome,stringsAsFactors=FALSE)
    fetch_genetic_markers$chromosome=chromosome
    
    reference_allele=paste0("\"", fetch_genetic_markers$reference_allele, "\"")
    reference_allele=as.data.frame(reference_allele,stringsAsFactors=FALSE)
    fetch_genetic_markers$reference_allele=reference_allele
    
    alternate_allele=paste0("\"", fetch_genetic_markers$alternate_allele, "\"")
    alternate_allele=as.data.frame(alternate_allele,stringsAsFactors=FALSE)
    fetch_genetic_markers$alternate_allele=alternate_allele
    
    marker_type=paste0("\"", fetch_genetic_markers$marker_type, "\"")
    marker_type=as.data.frame(marker_type,stringsAsFactors=FALSE)
    fetch_genetic_markers$marker_type=marker_type
    
    fetch_genetic_markers <- data.frame(marker_id=fetch_genetic_markers$marker_id,rs_id=fetch_genetic_markers$rs_id,chromosome=fetch_genetic_markers$chromosome,position_hg19=fetch_genetic_markers$position_hg19,position_hg38=fetch_genetic_markers$position_hg38,position_cM=fetch_genetic_markers$position_cM,reference_allele=fetch_genetic_markers$reference_allele,alternate_allele=fetch_genetic_markers$alternate_allele,marker_type=fetch_genetic_markers$marker_type,maf_gnomad_ALL=fetch_genetic_markers$maf_gnomad_ALL,maf_gnomad_AFR=fetch_genetic_markers$maf_gnomad_AFR,maf_gnomad_NFE=fetch_genetic_markers$maf_gnomad_NFE,maf_gnomad_FIN=fetch_genetic_markers$maf_gnomad_FIN,maf_gnomad_AMR=fetch_genetic_markers$maf_gnomad_AMR,maf_gnomad_EAS=fetch_genetic_markers$maf_gnomad_EAS,maf_gnomad_SAS=fetch_genetic_markers$maf_gnomad_SAS, stringsAsFactors=FALSE)
    current_genetic_markers=fetch_genetic_markers # current database markers
    
    #get markers common to current database
    genetic_markers$type=interaction(genetic_markers$chromosome,genetic_markers$position_hg19,genetic_markers$reference_allele,genetic_markers$alternate_allele)
    current_genetic_markers$type=interaction(current_genetic_markers$chromosome,current_genetic_markers$position_hg19,current_genetic_markers$reference_allele,current_genetic_markers$alternate_allele)
    common_markers=subset(genetic_markers,genetic_markers$type %in% current_genetic_markers$type)
    current_common_markers=subset(current_genetic_markers,current_genetic_markers$type %in% genetic_markers$type)
    
    common_markers_genotypes=subset(genotypes,genotypes$marker_id %in% common_markers$marker_id)
    common_markers_genotypes$marker_id=current_common_markers$marker_id
    
    common_markers_row_names=row.names(common_markers)
    
    
    new_markers=subset(genetic_markers,!row.names(genetic_markers) %in% common_markers_row_names)
    if(nrow(new_markers)+nrow(common_markers)!=nrow(genetic_markers)){stop("new markers are not labeled properly in marker_id column in table genetic_markers")}
    new_markers_genotypes=subset(genotypes,genotypes$marker_id %in% new_markers$marker_id)
    
    if(nrow(new_markers)>0)
    {
      new_markers$marker_id=(end_marker+1):(end_marker+nrow(new_markers))
      new_markers_genotypes$marker_id=(end_marker+1):(end_marker+nrow(new_markers))
    }
    
    ##check if individuals there
    Individuals=dbSendQuery(db, "SELECT * FROM Individuals;") # can not add LIMIT here as in SQL
    Individuals <- dbFetch(Individuals,)
    if(individuals$individual_id %in% Individuals$individual_id){
      stop("individual_id already exist in the database")
    }
    
    for (ii in seq_len(nrow(individuals))) {
      ii_command <- paste0("INSERT INTO Individuals (", paste(names(individuals), collapse=","), ") VALUES(", paste(individuals[ii, ], collapse=","), ");\n")
      cat(ii_command, file=mysql_commands, append=TRUE)
    }
    
    Samples=dbSendQuery(db, "SELECT * FROM Samples;") # can not add LIMIT here as in SQL
    Samples <- dbFetch(Samples,)
    
    if(samples$sample_id %in% Samples$sample_id  & individuals$individual_id %in% Individuals$individual_id){
      stop("sample_id already exist the database for the same individual_id")
    }
    
    for (ii in seq_len(nrow(samples))) {
      ii_command <- paste0("INSERT INTO Samples (", paste(names(samples), collapse=","), ") VALUES(", paste(samples[ii, ], collapse=","), ");\n")
      cat(ii_command, file=mysql_commands, append=TRUE)
    }
    
    ##check if pathogenic_mutations is there
    
    PathogenicMutations=dbSendQuery(db, "SELECT * FROM DiseaseCausingVariants;") # can not add LIMIT here as in SQL
    PathogenicMutations <- dbFetch(PathogenicMutations,)
    
    if(!(pathogenic_mutations$DCV_id %in% PathogenicMutations$mutation_id))
    {
      
      for (ii in seq_len(nrow(pathogenic_mutations))) {
        ii_command <- paste0("INSERT INTO DiseaseCausingVariants (", paste(names(pathogenic_mutations), collapse=","), ") VALUES(", paste(pathogenic_mutations[ii, ], collapse=","), ");\n")
        cat(ii_command, file=mysql_commands, append=TRUE)
      }
    }
    
    for (ii in seq_len(nrow(individuals_with_known_mutations))) {
      ii_command <- paste0("INSERT INTO IndividualsWithDiseaseCausingVariants (", paste(names(individuals_with_known_mutations), collapse=","), ") VALUES(", paste(individuals_with_known_mutations[ii, ], collapse=","), ");\n")
      cat(ii_command, file=mysql_commands, append=TRUE)
    }
    ##write ONLY common_markers_genotypes to sql file
    
    for (ii in seq_len(nrow(common_markers_genotypes))) {
      ii_command <- paste0("INSERT INTO Genotypes (", paste(names(common_markers_genotypes), collapse=","), ") VALUES(", paste(common_markers_genotypes[ii, ], collapse=","), ");\n")
      cat(ii_command, file=mysql_commands, append=TRUE)
    }
    
    
    ##write ONLY new_markers and new_markers_genotypes to the .sql file because common_markers_genotypes and common_markers are already in teh database
    
    if(nrow(new_markers)>0)
    {
      which(colnames(new_markers)=="type")
      new_markers=new_markers[,-which(colnames(new_markers)=="type")]
      
      for (ii in seq_len(nrow(new_markers))) {
        ii_command <- paste0("INSERT INTO GeneticMarkers (", paste(names(new_markers), collapse=","), ") VALUES(", paste(new_markers[ii, ], collapse=","), ");\n")
        cat(ii_command, file=mysql_commands, append=TRUE)
      }
      
      
      for (ii in seq_len(nrow(new_markers_genotypes))) {
        ii_command <- paste0("INSERT INTO Genotypes (", paste(names(new_markers_genotypes), collapse=","), ") VALUES(", paste(new_markers_genotypes[ii, ], collapse=","), ");\n")
        cat(ii_command, file=mysql_commands, append=TRUE)
      }
      
    }
  }
  # writing the sql commands to import second disease onwards into save_SQL_FILE
  
  print(paste0("SQL script is saved in ",save_SQL_FILE))
}
