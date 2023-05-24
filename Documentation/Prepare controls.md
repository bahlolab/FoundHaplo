1000 Genomes phase 3 haplotypes are used as controls. For every disease haplotype, control cohort should be run the same way a test cohort is run. 
Make sure you have downloaded and prepared 1000 Genomes phase 3 haplotypes and the sample IDs using [FoundHaplo/scripts/prepare_inputs/Initialise.sh](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/prepare_inputs/Initialise.sh).

Follow below steps everytime you test a new disease-causing variant using FoundHaplo.

1. Run [FoundHaplo/scripts/prepare_inputs/create_1000G_control_haplotypes.sh](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/prepare_inputs/create_1000G_control_haplotypes.sh) once for one disease-causing variant.
```bash
FoundHaplo_DIR=/mypath/FoundHaplo
DCV=DCV

$FoundHaplo_DIR/scripts/prepare_inputs/create_1000G_control_haplotypes.sh $FoundHaplo_DIR $DCV
```
Run example below,

```bash
$FoundHaplo_DIR/scripts/prepare_inputs/create_1000G_control_haplotypes.sh $FoundHaplo_DIR FAME1.chr8.119379052
```

Set the variables as below,,

* FoundHaplo_DIR : FoundHaplo directory i.e /mypath/FoundHaplo
* DCV : Name the disease variant of interest in the format of disease.chr.position. i.e FAME1.chr8.119379052. Use OMIM abbreviation for the disease.

The create_1000G_control_haplotypes.sh will create separate VCF files for each disease variant that you want to test inside five sub folders in FoundHaplo/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant corresponding to five super populations.   
Each sub folder must only contain VCF files with samples from the respective population listed in FoundHaplo/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/EUR.txt,AMR.txt,EAS.txt,SAS.txt and AFR.txt. 

Go back to the [documentation](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Guide%20to%20run%20FoundHaplo.md).
