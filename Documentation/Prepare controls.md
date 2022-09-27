1000 Genomes phase 3 haplotypes are used as controls. For every disease haplotype, control cohort should be run the same way a test cohort is run. 
Make sure you have downloaded 1000 Genomes phase 3 haplotypes and the sample IDs using FoundHaplo/scripts/prepare_inputs/Initialise.sh.

1. Run FoundHaplo/scripts/prepare_inputs/create_1000G_control_haplotypes.sh.
```bash
FoundHaplo/scripts/prepare_inputs/create_1000G_control_haplotypes.sh MAIN_PATH DCV CHROMOSOME VCFTOOLS_PATH
```

Arguments are explained below,

* MAIN_PATH : Path to FoundHaplo directory i.e path/FoundHaplo
* DCV : Name the disease variant of interest in the format of disease.chr.position. i.e FAME1.chr8.119379052.
* CHROMOSOME : Chromosome relevant to the interested disease variant without "chr" prefix
* VCFTOOLS_PATH : Path to vcftoools executable 

The create_1000G_control_haplotypes.sh will create separate VCF files for each disease variant that you want to test inside five sub folders corresponding to five super populations, each sub folder must only contain VCF files with samples from the respective population listed in FoundHaplo/input_files/public_data/1000G_haplotypes/EUR.txt,AMR.txt,EAS.txt,SAS.txt and AFR.txt. 

