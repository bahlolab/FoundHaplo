
Clone the repository
```bash
git clone --depth=1 https://github.com/bahlolab/Foundhaplo.git
```
Publicly available reference files required to run FoundHaplo

Use [FoundHaplo/scripts/prepare_inputs/Initialise.sh](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/prepare_inputs/Initialise.sh) as a guide to download the public data files into the relevant directories and software tools required by FoundHaplo

* [1000 Genomes phased 3 bed/bim/fam files](https://figshare.com/articles/dataset/1000_genomes_phase_3_files_with_SNPs_in_common_with_HapMap3/9208979?file=17838962) to harmonize input VCF data. 
* [1000 Genomes phase 3 haplotypes](http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/) are used as the control cohort when running FoundHaplo.
* [1000 Genomes phase 3 sample names by super population](https://www.internationalgenome.org/data-portal/data-collection/phase-3) to subset the control cohort to specifc super populations when running FoundHaplo.
* [Hapmap recombination files](https://github.com/bahlolab/FoundHaplo/tree/main/input_files/public_data/genetic_map_HapMapII_GRCh37) in hg19 to calculate the recombination rates.

Publicly available software tools required to run FoundHaplo

* [GenotypeHarmonizer](https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer-Download) tool is used to harmomize input VCF data to 1000Genomes.
* [Plink](https://zzz.bwh.harvard.edu/plink/plink2.shtml) version 1.9 is used to perform quality control steps on input VCF data. 
* [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/user-guide/download/) tool is used to annotate population frequencies form gnomAD data.
* [VCFtools](https://vcftools.github.io/downloads.html) version v0.1.13 and [BCFtools](http://www.htslib.org/download/) version 1.16 to perform queries on VCF files 
* [Nextflow](https://www.nextflow.io/) to parallely run FoundHaplo

Go back to the [main page](https://github.com/bahlolab/FoundHaplo).
