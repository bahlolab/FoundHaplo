# FoundHaplo

Tool to Identify individuals with inherited disease causing variants using SNP chip data.

Individuals who inherit the same genetic mutation from a common ancestor also share genomic regions either side of the shared disease-causing variant. This suggests that the presence of a disease-causing genetic variant can be inferred by assessing identity by descent sharing of the variant-associated haplotype between an individual known to have the disease-causing mutation and a patient with unknown aetiology. 

Using a Hidden Markov Model, Foundhaplo calculates the FH score, which is the likelihood of identity by descent between the samples of interest and a known disease causing haplotype and evaluate the resulted scores against a control cohort to predict individuals that carry disease causing haplotypes using SNP chip data. 

## A guide to use FoundHaplo

see [documentation](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Guide%20to%20run%20FoundHaplo.md) about input file formats and guide to use parallel processing.


## Install the R package using devtools packages

```R
devtools::install_github("bahlolab/FoundHaplo")
```

## Clone the repository
```bash
git clone --depth=1 https://github.com/bahlolab/FoundHaplo.git
```
