

Create_SQL_script_to_import(disease_hap_file,save_SQL_file,port,host,password,dbname,unix.socket,family_id,individual_id,father_id,mother_id,sex,ethnicity,ethnicity_superpopulation,ethnicity_method,sample_id,data_type,external_lab_id,impute_method,impute_panel,import_date,mutation_id,disease,disease_id,omim_id,gene=inheritance_model,chr,start_position_hg19,end_position_hg19,start_position_hg38,end_position_hg38,start_position_cM,end_position_cM,genotype,validated,validation_method,validation_note)
{
  
  library(data.table)
  library(DBI)
  library(RMariaDB)
  disease_hap_file <-fread(disease_hap_file,skip = "#CHROM")
  
  options(scipen=99)
  
  individuals <- data.frame(family_id=family_id, individual_id=individual_id, father_id=father_id, mother_id=mother_id,
                            sex=sex, ethnicity=paste0("\"",ethnicity,"\""),
                            ethnicity_superpopulation=paste0("\"",ethnicity_superpopulation,"\""),ethnicity_method=paste0("\"",ethnicity_method,"\""),
                            stringsAsFactors=FALSE)
  
  
  samples <- data.frame(sample_id=sample_id, individual_id=individual_id,data_type=paste0("\"",data_type,"\""),external_lab_id=paste0("\"",external_lab_id,"\""),impute_method=paste0("\"",impute_method,"\""),impute_panel=paste0("\"",impute_panel,"\""),import_date=paste0("\"",import_date,"\""),
                        stringsAsFactors=FALSE)
  
  pathogenic_mutations <- data.frame(mutation_id=mutation_id,
                                     disease=paste0("\"",disease,"\""), disease_id=paste0("\"",disease_id,"\""),
                                     omim_id=omim_id, gene=paste0("\"",gene,"\""), inheritance_model=paste0("\"",inheritance_model,"\""),
                                     chr=paste0("\"",chr,"\""), start_position_hg19=start_position_hg19, end_position_hg19=end_position_hg19,start_position_hg38=start_position_hg38,end_position_hg38=end_position_hg38,
                                     start_position_cM=start_position_cM, end_position_cM=end_position_cM,
                                     stringsAsFactors=FALSE)
  
  individuals_with_known_mutations <- data.frame(individual_id=individual_id, mutation_id=mutation_id, genotype=genotype,
                                                 validated=validated, validation_method=paste0("\"",validation_method,"\""),
                                                 validation_note=paste0("\"",validation_note,"\""),
                                                 stringsAsFactors=FALSE)
  
  fields=c("R2","AF","AF_raw", "AF_afr", "AF_sas","AF_amr", "AF_eas", "AF_nfe", "AF_fin")
  
  get_info_fields <- function(info, fields) {
    info_split <- strsplit(strsplit(info, ";")[[1]], "=")
    info_split_names <- sapply(info_split, "[", 1)
    info_split_names[1:6]=c("R2","AF","AF_raw", "AF_afr", "AF_sas","AF_amr", "AF_eas", "AF_nfe", "AF_fin")
    info_split_values <- sapply(info_split, "[", 2)
    
    return(info_split_values[match(fields, fields)])
  }
  
  hap_info_extract <- do.call(rbind, lapply(disease_hap_file$INFO, function(x) get_info_fields(x, fields)))
  hap_info_extract <- apply(hap_info_extract, 2, as.numeric)
  hap_info_extract <- apply(hap_info_extract, 2, function(x) ifelse(is.na(x), 0, x))
  hap_info_extract <- as.data.frame(hap_info_extract, stringsAsFactors=FALSE)
  rownames(hap_info_extract) <- rownames(disease_hap_file)
  colnames(hap_info_extract) <- fields
  
  disease_hap_file=as.data.frame(disease_hap_file)
  
  genetic_markers <- data.frame(marker_id=1:nrow(disease_hap_file), rs_id=paste0("\"", disease_hap_file$ID, "\""),
                                chr=paste0("\"chr", disease_hap_file[,"#CHROM"], "\""),
                                position_hg19=disease_hap_file$POS,
                                position_hg38=-99999,
                                position_cM=-99999,
                                ref=paste0("\"", disease_hap_file$REF, "\""),
                                alt=paste0("\"", disease_hap_file$ALT, "\""),
                                marker_type="\"SNP\"",
                                maf_gnomad_ALL=as.numeric(hap_info_extract$AF_raw),
                                maf_gnomad_AFR=as.numeric(hap_info_extract$AF_afr),
                                maf_gnomad_NFE=as.numeric(hap_info_extract$AF_nfe),
                                maf_gnomad_FIN=as.numeric(hap_info_extract$AF_fin),
                                maf_gnomad_AMR=as.numeric(hap_info_extract$AF_amr),
                                maf_gnomad_EAS=as.numeric(hap_info_extract$AF_eas),
                                maf_gnomad_SAS=as.numeric(hap_info_extract$AF_sas),
                                stringsAsFactors=FALSE)
  
  
  R2=sapply(strsplit(disease_hap_file$INFO,";",fixed=TRUE),"[[", 1)
  table(R2 %like% "R2")
  
  R2=sapply(strsplit(R2,"=",fixed=TRUE),"[[", 2)
  
  
  
  genotypes <- data.frame(marker_id=1:nrow(disease_hap_file), sample_id=sample_id, genotype=disease_hap_file$h1,imputed=1,imputation_quality=R2,
                          stringsAsFactors=FALSE)
  
  genotypes <- genotypes[!is.na(genotypes$genotype), ]
  dummy = strsplit(genotypes$genotype,":")[[1]]
  
  dummy=sapply(strsplit(genotypes$genotype,":",fixed=TRUE),"[[", 1)
  dummy=sapply(strsplit(dummy,"/",fixed=TRUE),"[[", 1)
  
  genotypes$genotype=dummy
  
  mysql_commands <- save_SQL_file
  cat("USE FoundHaploDB;\n", file=mysql_commands)
  
  db = dbConnect(RMariaDB::MariaDB(),bigint = 'integer',port=port,host=host,user ='remote_usr',password=password,dbname=dbname,unix.socket=unix.socket)
  
  fetch_genetic_markers=dbSendQuery(db, "SELECT * FROM GeneticMarkers;") # can not add LIMIT here as in SQL
  fetch_genetic_markers <- dbFetch(fetch_genetic_markers,)
  
  end_marker=max(fetch_genetic_markers$marker_id)
  
  rs_id=paste0("\"", fetch_genetic_markers$rs_id, "\"")
  rs_id=as.data.frame(rs_id,stringsAsFactors=FALSE)
  
  fetch_genetic_markers$rs_id=rs_id
  
  chr=paste0("\"", fetch_genetic_markers$chr, "\"")
  chr=as.data.frame(chr,stringsAsFactors=FALSE)
  fetch_genetic_markers$chr=chr
  
  ref=paste0("\"", fetch_genetic_markers$ref, "\"")
  ref=as.data.frame(ref,stringsAsFactors=FALSE)
  fetch_genetic_markers$ref=ref
  
  alt=paste0("\"", fetch_genetic_markers$alt, "\"")
  alt=as.data.frame(alt,stringsAsFactors=FALSE)
  fetch_genetic_markers$alt=alt
  
  marker_type=paste0("\"", fetch_genetic_markers$marker_type, "\"")
  marker_type=as.data.frame(marker_type,stringsAsFactors=FALSE)
  fetch_genetic_markers$marker_type=marker_type
  
  fetch_genetic_markers <- data.frame(marker_id=fetch_genetic_markers$marker_id,rs_id=fetch_genetic_markers$rs_id,chr=fetch_genetic_markers$chr,position_hg19=fetch_genetic_markers$position_hg19,position_hg38=fetch_genetic_markers$position_hg38,position_cm=fetch_genetic_markers$position_cm,ref=fetch_genetic_markers$ref,alt=fetch_genetic_markers$alt,marker_type=fetch_genetic_markers$marker_type,maf_gnomad_ALL=fetch_genetic_markers$maf_gnomad_ALL,maf_gnomad_AFR=fetch_genetic_markers$maf_gnomad_AFR,maf_gnomad_NFE=fetch_genetic_markers$maf_gnomad_NFE,maf_gnomad_FIN=fetch_genetic_markers$maf_gnomad_FIN,maf_gnomad_AMR=fetch_genetic_markers$maf_gnomad_AMR,maf_gnomad_EAS=fetch_genetic_markers$maf_gnomad_EAS,maf_gnomad_SAS=fetch_genetic_markers$maf_gnomad_SAS, stringsAsFactors=FALSE)
  
  
  
  current_genetic_markers=fetch_genetic_markers # current database markers
  
  
  #run the import script for the new file to import
  #get markers common to current database
  genetic_markers$type=interaction(genetic_markers$chr,genetic_markers$position_hg19,genetic_markers$ref,genetic_markers$alt)
  current_genetic_markers$type=interaction(current_genetic_markers$chr,current_genetic_markers$position_hg19,current_genetic_markers$ref,current_genetic_markers$alt)
  common_markers=subset(genetic_markers,genetic_markers$type %in% current_genetic_markers$type)
  current_common_markers=subset(current_genetic_markers,current_genetic_markers$type %in% genetic_markers$type)
  
  common_markers_genotypes=subset(genotypes,genotypes$marker_id %in% common_markers$marker_id)
  common_markers_genotypes$marker_id=current_common_markers$marker_id
  #write ONLY common_markers_genotypes to sql file
  
  #get new markers current database
  common_markers_row_names=row.names(common_markers)
  
  
  new_markers=subset(genetic_markers,!row.names(genetic_markers) %in% common_markers_row_names)
  nrow(new_markers)+nrow(common_markers)==nrow(genetic_markers)
  new_markers_genotypes=subset(genotypes,genotypes$marker_id %in% new_markers$marker_id)
  
  if(nrow(new_markers)>0)
  {
    new_markers$marker_id=(end_marker+1):(end_marker+nrow(new_markers))
    new_markers_genotypes$marker_id=(end_marker+1):(end_marker+nrow(new_markers))
  }
  #### run carefully
  mysql_commands <- save_SQL_file
  cat("USE FoundHaploDB;\n", file=mysql_commands)
  
  ##check if individuals there
  for (ii in seq_len(nrow(individuals))) {
    ii_command <- paste0("INSERT INTO Individuals (", paste(names(individuals), collapse=","), ") VALUES(", paste(individuals[ii, ], collapse=","), ");\n")
    cat(ii_command, file=mysql_commands, append=TRUE)
  }
  
  for (ii in seq_len(nrow(samples))) {
    ii_command <- paste0("INSERT INTO Samples (", paste(names(samples), collapse=","), ") VALUES(", paste(samples[ii, ], collapse=","), ");\n")
    cat(ii_command, file=mysql_commands, append=TRUE)
  }
  
  ##check if pathogenic_mutations is there
  
  for (ii in seq_len(nrow(pathogenic_mutations))) {
    ii_command <- paste0("INSERT INTO PathogenicMutations (", paste(names(pathogenic_mutations), collapse=","), ") VALUES(", paste(pathogenic_mutations[ii, ], collapse=","), ");\n")
    cat(ii_command, file=mysql_commands, append=TRUE)
  }
  
  for (ii in seq_len(nrow(individuals_with_known_mutations))) {
    ii_command <- paste0("INSERT INTO IndividualsWithKnownMutations (", paste(names(individuals_with_known_mutations), collapse=","), ") VALUES(", paste(individuals_with_known_mutations[ii, ], collapse=","), ");\n")
    cat(ii_command, file=mysql_commands, append=TRUE)
  }
  ##write ONLY common_markers_genotypes to sql file
  
  for (ii in seq_len(nrow(common_markers_genotypes))) {
    ii_command <- paste0("INSERT INTO Genotypes (", paste(names(common_markers_genotypes), collapse=","), ") VALUES(", paste(common_markers_genotypes[ii, ], collapse=","), ");\n")
    cat(ii_command, file=mysql_commands, append=TRUE)
  }
  
  
  ##write ONLY new_markers and new_markers_genotypes to sql file
  
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
  print("SQL script is saved in ",save_SQL_file) 
}
