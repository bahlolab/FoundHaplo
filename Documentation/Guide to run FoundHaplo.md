1. [Overview of FoundHaplo](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Overview%20of%20FoundHaplo.md)

### Prepare input data

3. [Set-up software tools and publicly available input files](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Publicly%20available%20Input%20files%20and%20software%20tools.md)
4. [Prepare known disease haplotypes](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Prepare%20known%20disease%20haplotypes.md) 
     * [Prepare controls](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Prepare%20controls.md)
5. [Prepare test cohorts](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Prepare%20test%20samples.md)

### Run FoundHaplo

Details of all the functions in FoundHaplo R package can be found [here](https://github.com/bahlolab/FoundHaplo/blob/main/vignettes).

7. [Parameters that must be specified by the user to generate the IBD report for a given disease-causing variant](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Parameters%20in%20the%20algorithm.md)
8. [Run FoundHaplo parallely using Nextflow pipeline](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Parallel%20processing.md)

### Analyse FoundHaplo results

9. Run [Analyse_FH.R](https://github.com/bahlolab/FoundHaplo/blob/main/R/Analyse_FH.R) to predict test samples likely to carry the tested disease-causing variants.

```R
Analyse_FH(path_results="FoundHaplo/results/FH_IBD_scores/results.txt",path_to_save_FH_output="FoundHaplo/results/FH_Analysis",critical_percentile=0.99)
```

All the parameters that user has to specify are described below

1. **path_results** Path to a single .txt file with all the FH scores
2. **path_to_save_FH_output** Path to save the graphical output of the FH scores 
3. **critical_percentile** Critical percentile of the control cohort to derive predictions. Recommend above 0.999 for large cohorts like UKBB

Go back to the [documentation](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Guide%20to%20run%20FoundHaplo.md).
