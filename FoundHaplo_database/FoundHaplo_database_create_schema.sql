CREATE DATABASE FoundHaploDB;

USE FoundHaploDB;

CREATE TABLE `Individuals` (
    `family_id` int NOT NULL,
    `individual_id` int NOT NULL AUTO_INCREMENT,
    `father_id` int,
    `mother_id` int,
    `sex` tinyint,
    `sex_method` varchar(255),
    `ancestral_population` varchar(255),
    `ancestral_superpopulation` varchar(255) NOT NULL,
    `ancestry_method` varchar(255) NOT NULL,
    PRIMARY KEY (`individual_id`)
);

CREATE TABLE `Samples` (
    `sample_id` int NOT NULL AUTO_INCREMENT,
    `individual_id` int NOT NULL,
    `data_type` varchar(255) NOT NULL,
    `external_lab_id` varchar(255) NOT NULL,
    `external_source` varchar(255),
    `phasing_method` varchar(255) NOT NULL,
    `impute_method` varchar(255) NOT NULL,
    `impute_phasing_panel` varchar(255) NOT NULL,
    `import_date` date NOT NULL,
    PRIMARY KEY (`sample_id`),
    FOREIGN KEY (`individual_id`) 
        REFERENCES `Individuals`(`individual_id`)
);

CREATE TABLE `DiseaseCausingVariants` (
    `DCV_id` int NOT NULL AUTO_INCREMENT,
    `disease_name` varchar(255) NOT NULL,
    `omim_id` int,
    `gene` varchar(255),
    `genomic_region` varchar(255) NOT NULL,
    `inheritance_model` varchar(255) NOT NULL,
    `chromosome` varchar(15) NOT NULL,
    `start_position_hg19` int NOT NULL,
    `end_position_hg19` int,
    `start_position_hg38` int,
    `end_position_hg38` int,
    `start_position_cM` double,
    `end_position_cM` double,
    PRIMARY KEY (`mutation_id`)
);

CREATE TABLE `IndividualsWithDiseaseCausingVariants` (
    `individual_id` int NOT NULL,
    `DCV_id` int NOT NULL,
    `genotype` tinyint NOT NULL,
    `validated` BOOLEAN NOT NULL,
    `validation_method` varchar(255) NOT NULL,
    `validation_note` varchar(255),
    FOREIGN KEY (`individual_id`) 
        REFERENCES `Individuals`(`individual_id`),
    FOREIGN KEY (`DCV_id`) 
        REFERENCES `DiseaseCausingVariants`(`DCV_id`)
);

CREATE TABLE `GeneticMarkers` (
    `marker_id` bigint NOT NULL AUTO_INCREMENT,
    `rs_id` varchar(31),
    `chromosome` varchar(15) NOT NULL,
    `position_hg19` int NOT NULL,
    `position_hg38` int,
    `position_cM` double,
    `reference_allele` varchar(255) NOT NULL,
    `alternate_allele` varchar(255) NOT NULL,
    `marker_type` varchar(255) NOT NULL,
    `maf_gnomad_ALL` FLOAT NOT NULL,
    `maf_gnomad_AFR` FLOAT NOT NULL,
    `maf_gnomad_NFE` FLOAT NOT NULL,
    `maf_gnomad_FIN` FLOAT NOT NULL,
    `maf_gnomad_AMR` FLOAT NOT NULL,
    `maf_gnomad_EAS` FLOAT NOT NULL,
    `maf_gnomad_SAS` FLOAT NOT NULL,
    PRIMARY KEY (`marker_id`)
);

CREATE TABLE `Genotypes` (
    `marker_id` bigint NOT NULL,
    `sample_id` int NOT NULL,
    `genotype` tinyint NOT NULL,
    `imputed` BOOLEAN NOT NULL,
    `imputation_quality` FLOAT NOT NULL,
    PRIMARY KEY (`marker_id`, `sample_id`), 
    FOREIGN KEY (`marker_id`) 
        REFERENCES `GeneticMarkers`(`marker_id`),
    FOREIGN KEY (`sample_id`) 
        REFERENCES `Samples`(`sample_id`)
);



