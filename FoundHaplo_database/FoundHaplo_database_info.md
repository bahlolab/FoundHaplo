FoundHaplo_database_info
================
Erandee
2023-06-06

FoundHaplo database schema consists of six relational tables as listed
below. Tables are automatically created by the script in
FoundHaplo/FoundHaplo_database/FoundHaplo_database_create_schema.sql

## R Markdown

|        Field name         | Data type |                               Description                                |   Key   |
|:-------------------------:|:---------:|:------------------------------------------------------------------------:|:-------:|
|         family_id         |    int    |   Unique family identifier (One individual can have multiple samples)    |         |
|       individual_id       |    int    |                       Unique individual identifier                       | Primary |
|         father_id         |    int    |                   Individual identifier of the father                    |         |
|         mother_id         |    int    |                   Individual identifier of the mother                    |         |
|            sex            |  tinyint  |                             Male=1, Female=2                             |         |
|        sex_method         |  varchar  |                       Method of estimating the sex                       |         |
|   ancestral_population    |  varchar  |    Inferred ancestral population based on 1000 Genomes classification    |         |
| ancestral_superpopulation |  varchar  | Inferred ancestral super-population based on 1000 Genomes classification |         |
|      ancestry_method      |  varchar  |             Method of inferring ancestry ex: PCA or reported             |         |

Individuals table

|      Field name      | Data type |                             Description                             |              Key               |
|:--------------------:|:---------:|:-------------------------------------------------------------------:|:------------------------------:|
|      sample_id       |    int    | Unique sample identifier (One individual can have multiple samples) |            Primary             |
|    individual_id     |    int    |                    Unique individual identifier                     | Foreign from table Individuals |
|      data_type       |  varchar  |              Type of data ex: WGS/SNP genotyping etc.               |                                |
|   external_lab_id    |  varchar  |                  Unique original sample identifier                  |                                |
|   external_source    |  varchar  |                     Source of sample collection                     |                                |
|    phasing_method    |  varchar  |          Method used for phasing ex: trio, duo or related           |                                |
|    impute_method     |  varchar  |             Method used for imputation ex: MIS, TOPMed              |                                |
| impute_phasing_panel |  varchar  |         Reference panel used for imputing ex: 1000 Genomes          |                                |
|     import_date      |  varchar  |          The date the sample was imported to the database           |                                |

Samples table

|     Field name      | Data type |                                      Description                                       |              Key               |
|:-------------------:|:---------:|:--------------------------------------------------------------------------------------:|:------------------------------:|
|       DCV_id        |    int    |                       Unique identifier for a specific mutation                        |            Primary             |
|    disease_name     |  varchar  | Unique OMIM abbreviation (This should be the same as the first part of DCV parameter)  | Foreign from table Individuals |
|       omim_id       |    int    |                Unique OMIM identifier for the disease variant/mutation                 |                                |
|        gene         |  varchar  |                              Gene containing the mutatio                               |                                |
|   genomic_region    |  varchar  |                                 Exonic, intronic etc.                                  |                                |
|  inheritance_model  |  varchar  | Autosomal dominant (AD), recessive (AR) or X-linked dominant or recessive (XLD or XLR) |                                |
|     chromosome      |  varchar  |                                  Chromosome ex: chr1                                   |                                |
| start_position_hg19 |    int    |                    Start base pair position of the variant in hg19                     |                                |
|  end_position_hg19  |    int    |                     End base pair position of the variant in hg19                      |                                |
| start_position_hg38 |    int    |                    EStart base pair position of the variant in hg38                    |                                |
|  end_position_hg38  |    int    |                      End base pair position of the variant in hg3                      |                                |
|  start_position_cM  |  double   |                          Start position of the variant in cM                           |                                |
|   end_position_cM   |  double   |                           End position of the variant in cM                            |                                |

DiseaseCausingVariants table

|    Field name     | Data type |                          Description                           |                    Key                    |
|:-----------------:|:---------:|:--------------------------------------------------------------:|:-----------------------------------------:|
|   individual_id   |    int    |           Unique identifier for a specific mutation            |      Foreign from table Individuals       |
|      DCV_id       |    int    |           Unique identifier for a specific mutation            | Foreign from table DiseaseCausingVariants |
|     genotype      |    int    | INumber of copies of the identified disease variant ex: 1 or 2 |                                           |
|     validated     |  BOOLEAN  |               Is the disease variant validated?                |                                           |
| validation_method |  varchar  |        Validation method ex: RP-PCR, bioinformatic etc.        |                                           |
|  validation_note  |  varchar  |                Extra note on validation, if any                |                                           |

IndividualsWithDiseaseCausingVariants table

|    Field name    | Data type |                     Description                     |   Key   |
|:----------------:|:---------:|:---------------------------------------------------:|:-------:|
|    marker_id     |  bigint   |        Unique identifier of genetic markers         | Primary |
|      rs_id       |  varchar  |              Reference SNP identifier               |         |
|    chromosome    |    int    |                 Chromosome ex: chr1                 |         |
|  position_hg19   |    int    |             Base pair position in hg19              |         |
|  position_hg38   |    int    |             Base pair position in hg38              |         |
|   position_cM    |  double   |                   position in cM                    |         |
| reference_allele |  varchar  |                  Reference allele                   |         |
| alternate_allele |  varchar  |                  Alternate allele                   |         |
|   marker_type    |  varchar  |      SNP (FoundHaplo currently uses only SNPs)      |         |
|  maf_gnomad_ALL  |   FLOAT   |              General allele frequency               |         |
|  maf_gnomad_AFR  |   FLOAT   |       Allele frequency in African population        |         |
|  maf_gnomad_NFE  |   FLOAT   | Allele frequency in Non-finnish European population |         |
|  maf_gnomad_FIN  |   FLOAT   |       Allele frequency in Finnish population        |         |
|  maf_gnomad_AMR  |   FLOAT   |       Allele frequency in American population       |         |
|  maf_gnomad_EAS  |   FLOAT   |      Allele frequency in East Asian population      |         |
|  maf_gnomad_SAS  |   FLOAT   |     Allele frequency in South Asian population      |         |

GeneticMarkers table

|     Field name     | Data type |                              Description                              |            Key             |
|:------------------:|:---------:|:---------------------------------------------------------------------:|:--------------------------:|
|     marker_id      |  bigint   |                 Unique identifier of genetic markers                  |          Primary           |
|     sample_id      |    int    |  Unique sample identifier (One individual can have multiple samples)  | Foreign from table Samples |
|      genotype      |  tinyint  | “0” denotes the reference allele and “1” denotes the alternate allele |                            |
|      imputed       |  BOOLEAN  |                        Is the marker imputed?                         |                            |
| imputation_quality |   FLOAT   |                       Imputation quality score                        |                            |

Genotypes table
