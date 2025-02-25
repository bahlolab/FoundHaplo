FoundHaplo database information
================
Erandee
2023-06-18

FoundHaplo database schema consists of six relational tables as listed
below. Tables are automatically created by the script in
FoundHaplo/FoundHaplo_database/FoundHaplo_database_create_schema.sql

1.  “data_type” in Samples table : Disease haplotypes should be SNP
    genotyped since FoundHaplo currently only works on SNP array data

2.  “phasing_method” in Samples table : Disease haplotypes should be
    phased using duo, trio or related samples so that disease haplotypes
    can be created by pedigree phasing

3.  “marker_type” in GeneticMarkers table : Genetic markers should be
    bi-allelic SNPs. FoundHaplo algorithm only uses bi-allelic SNPs for
    inferrence

4.  “ancestral_population” and “ancestral_superpopulation” in
    Individuals table : Ancestry should be classified based on 1000
    Genomes ancestry categories since 1000 Genomes data is used as the
    control cohort for inference

|        Field name         | Data type |                               Description                                |   Key   |
|:-------------------------:|:---------:|:------------------------------------------------------------------------:|:-------:|
|         family_id         |    int    |   Unique family identifier (One individual can have multiple samples)    |         |
|       individual_id       |    int    |                       Unique individual identifier                       | Primary |
|         father_id         |    int    |                   Individual identifier of the father                    |         |
|         mother_id         |    int    |                   Individual identifier of the mother                    |         |
|            sex            |  tinyint  |                             Male=1, Female=2                             |         |
|        sex_method         |  varchar  |                       Method of inferring the sex                        |         |
|   ancestral_population    |  varchar  |    Inferred ancestral population based on 1000 Genomes classification    |         |
| ancestral_superpopulation |  varchar  | Inferred ancestral super-population based on 1000 Genomes classification |         |
|      ancestry_method      |  varchar  |             Method of inferring ancestry eg: PCA or reported             |         |

Individuals table

|      Field name      | Data type |                                                   Description                                                   |              Key               |
|:--------------------:|:---------:|:---------------------------------------------------------------------------------------------------------------:|:------------------------------:|
|      sample_id       |    int    |                       Unique sample identifier (One individual can have multiple samples)                       |            Primary             |
|    individual_id     |    int    |                                          Unique individual identifier                                           | Foreign from table Individuals |
|      data_type       |  varchar  |                                   Type of data such as WGS/WES/SNP genotyping                                   |                                |
|   external_lab_id    |  varchar  |                                        Unique original sample identifier                                        |                                |
|   external_source    |  varchar  |                                           Source of sample collection                                           |                                |
|    phasing_method    |  varchar  | Method used for phasing such as trio, duo, related or using imputation tools such as Michigan Imputation Server |                                |
|    impute_method     |  varchar  |                                   Method used for imputation eg: MIS, TOPMed                                    |                                |
| impute_phasing_panel |  varchar  |                               Reference panel used for imputing eg: 1000 Genomes                                |                                |
|     import_date      |   date    |                                The date the sample was imported to the database                                 |                                |

Samples table

|     Field name      | Data type |                                      Description                                       |              Key               |
|:-------------------:|:---------:|:--------------------------------------------------------------------------------------:|:------------------------------:|
|       DCV_id        |    int    |                       Unique identifier for a specific mutation                        |            Primary             |
|    disease_name     |  varchar  | Unique OMIM abbreviation (This should be the same as the first part of DCV parameter)  | Foreign from table Individuals |
|       omim_id       |    int    |                     Unique OMIM identifier for the disease variant                     |                                |
|        gene         |  varchar  |                              Gene containing the mutation                              |                                |
|   genomic_region    |  varchar  |                                 Exonic, intronic etc.                                  |                                |
|  inheritance_model  |  varchar  | Autosomal dominant (AD), recessive (AR) or X-linked dominant or recessive (XLD or XLR) |                                |
|     chromosome      |  varchar  |                                  Chromosome eg: chr1                                   |                                |
| start_position_hg19 |    int    |                    Start base pair position of the variant in hg19                     |                                |
|  end_position_hg19  |    int    |                     End base pair position of the variant in hg19                      |                                |
| start_position_hg38 |    int    |                    EStart base pair position of the variant in hg38                    |                                |
|  end_position_hg38  |    int    |                      End base pair position of the variant in hg3                      |                                |
|  start_position_cM  |  double   |                          Start position of the variant in cM                           |                                |
|   end_position_cM   |  double   |                           End position of the variant in cM                            |                                |

DiseaseCausingVariants table

|    Field name     | Data type |                          Description                          |                          Key                          |
|:-----------------:|:---------:|:-------------------------------------------------------------:|:-----------------------------------------------------:|
|   individual_id   |    int    |                 Unique individual identifier                  |      Primary and Foreign from table Individuals       |
|      DCV_id       |    int    |           Unique identifier for a specific mutation           | Primary and Foreign from table DiseaseCausingVariants |
|     genotype      |    int    | Number of copies of the identified disease variant eg: 1 or 2 |                                                       |
|     validated     |  BOOLEAN  |               Is the disease variant validated?               |                                                       |
| validation_method |  varchar  |       Validation method eg: RP-PCR, bioinformatic etc.        |                                                       |
|  validation_note  |  varchar  |               Extra note on validation, if any                |                                                       |

IndividualsWithDiseaseCausingVariants table

|    Field name    | Data type |                            Description                            |   Key   |
|:----------------:|:---------:|:-----------------------------------------------------------------:|:-------:|
|    marker_id     |  bigint   |               Unique identifier of genetic markers                | Primary |
|      rs_id       |  varchar  |                     Reference SNP identifier                      |         |
|    chromosome    |    int    |                        Chromosome eg: chr1                        |         |
|  position_hg19   |    int    |                    Base pair position in hg19                     |         |
|  position_hg38   |    int    |                    Base pair position in hg38                     |         |
|   position_cM    |  double   |                          position in cM                           |         |
| reference_allele |  varchar  |                         Reference allele                          |         |
| alternate_allele |  varchar  |                         Alternate allele                          |         |
|   marker_type    |  varchar  | SNP/insertion/deletion etc. (FoundHaplo currently uses only SNPs) |         |
|  maf_gnomad_ALL  |   FLOAT   |                     General allele frequency                      |         |
|  maf_gnomad_AFR  |   FLOAT   |              Allele frequency in African population               |         |
|  maf_gnomad_NFE  |   FLOAT   |        Allele frequency in Non-finnish European population        |         |
|  maf_gnomad_FIN  |   FLOAT   |              Allele frequency in Finnish population               |         |
|  maf_gnomad_AMR  |   FLOAT   |              Allele frequency in American population              |         |
|  maf_gnomad_EAS  |   FLOAT   |             Allele frequency in East Asian population             |         |
|  maf_gnomad_SAS  |   FLOAT   |            Allele frequency in South Asian population             |         |

GeneticMarkers table

|     Field name     | Data type |                              Description                              |                      Key                      |
|:------------------:|:---------:|:---------------------------------------------------------------------:|:---------------------------------------------:|
|     marker_id      |  bigint   |                 Unique identifier of genetic markers                  | Primary and Foreign from table GeneticMarkers |
|     sample_id      |    int    |  Unique sample identifier (One individual can have multiple samples)  |    Primary and Foreign from table Samples     |
|      genotype      |  tinyint  | “0” denotes the reference allele and “1” denotes the alternate allele |                                               |
|      imputed       |  BOOLEAN  |                        Is the marker imputed?                         |                                               |
| imputation_quality |   FLOAT   |                       Imputation quality score                        |                                               |

Genotypes table
