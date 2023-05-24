CREATE DATABASE FoundHaploDB;

USE FoundHaploDB;

CREATE TABLE `Individuals` (
    `family_id` int NOT NULL,
    `individual_id` int NOT NULL AUTO_INCREMENT,
    `father_id` int,
    `mother_id` int,
    `sex` tinyint,
    `ethnicity` varchar(255),
    `ethnicity_superpopulation` varchar(255),
    `ethnicity_method` varchar(255),
    PRIMARY KEY (`individual_id`)
);

CREATE TABLE `Samples` (
    `sample_id` int NOT NULL AUTO_INCREMENT,
    `individual_id` int NOT NULL,
    `data_type` varchar(255) NOT NULL,
    `external_lab_id` varchar(255) NOT NULL,
    `impute_method` varchar(255) NOT NULL,
    `impute_panel` varchar(255),
    `import_date` varchar(255) NOT NULL,
    PRIMARY KEY (`sample_id`),
    FOREIGN KEY (`individual_id`) 
        REFERENCES `Individuals`(`individual_id`)
);

CREATE TABLE `PathogenicMutations` (
    `mutation_id` int NOT NULL AUTO_INCREMENT,
    `disease` varchar(255) NOT NULL,
    `disease_id` varchar(20) NOT NULL,
    `omim_id` int,
    `gene` varchar(255),
    `inheritance_model` varchar(255) NOT NULL,
    `chr` varchar(15) NOT NULL,
    `start_position_hg19` int NOT NULL,
    `end_position_hg19` int,
    `start_position_hg38` int,
    `end_position_hg38` int,
    `start_position_cM` double,
    `end_position_cM` double,
    PRIMARY KEY (`mutation_id`)
);

CREATE TABLE `IndividualsWithKnownMutations` (
    `individual_id` int NOT NULL,
    `mutation_id` int NOT NULL,
    `genotype` tinyint NOT NULL,
    `validated` BOOLEAN,
    `validation_method` varchar(255),
    `validation_note` varchar(255),
    FOREIGN KEY (`individual_id`) 
        REFERENCES `Individuals`(`individual_id`),
    FOREIGN KEY (`mutation_id`) 
        REFERENCES `PathogenicMutations`(`mutation_id`)
);

CREATE TABLE `GeneticMarkers` (
    `marker_id` bigint NOT NULL AUTO_INCREMENT,
    `rs_id` varchar(31),
    `chr` varchar(15) NOT NULL,
    `position_hg19` int NOT NULL,
    `position_hg38` int,
    `position_cM` double,
    `ref` varchar(255) NOT NULL,
    `alt` varchar(255) NOT NULL,
    `marker_type` varchar(255),
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
    `imputed` BOOLEAN,
    `imputation_quality` FLOAT,
    PRIMARY KEY (`marker_id`, `sample_id`), 
    FOREIGN KEY (`marker_id`) 
        REFERENCES `GeneticMarkers`(`marker_id`),
    FOREIGN KEY (`sample_id`) 
        REFERENCES `Samples`(`sample_id`)
);



