1. [Overview of FoundHaplo](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Overview%20of%20FoundHaplo.md)

### Prepare input data

2. [Set-up software tools and publicly available input files](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Publicly%20available%20Input%20files%20and%20software%20tools.md)
3. [Prepare known disease haplotypes](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Prepare%20known%20disease%20haplotypes.md) 
     * [Prepare controls](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Prepare%20controls.md)
4. [Prepare test cohorts](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Prepare%20test%20samples.md)

### Run FoundHaplo

Details of all the functions in FoundHaplo R package can be found [here](https://github.com/bahlolab/FoundHaplo/blob/main/vignettes).

5. Run FoundHaplo [parallely using Nextflow pipeline](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Parallel%20processing%20with%20Nextflow.md)

If you do not have the Nextflow pipeline, run the main R function [Generate_FH_score](https://github.com/bahlolab/FoundHaplo/blob/main/R/Generate_FH_score.R) sequentially as explained [here](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Parameters%20in%20the%20Generate_FH_score.md) 

### Analyse FoundHaplo results

6. After running FoundHaplo, Run [Analyse_FH.R](https://github.com/bahlolab/FoundHaplo/blob/main/R/Analyse_FH.R) to get the predictions and plot the test samples likely to carry the tested disease-causing variants.

```R
Analyse_FH(results_FILE="/mypath/FoundHaplo/results/FH_IBD_scores/results.txt",save_FH_output_DIR="/mypath/FoundHaplo/results/FH_Analysis",critical_percentile=99,from_control = FALSE)
```

All the parameters that the user has to specify are described below

* results_FILE = Path to a single .txt file with all the FH scores
* save_FH_output_DIR = Directory to save the graphical output of the FH scores 
* critical_percentile = Critical percentile of the control cohort (or the test cohort itself) to derive predictions. Recommend above the 90th percentile; default is 99th percentile
* from_control = TRUE if the critical value is calculated from a control cohort or FALSE otherwise. Recommend running a control cohort; default is, from_control = TRUE

Go back to the [main page](https://github.com/bahlolab/FoundHaplo).
