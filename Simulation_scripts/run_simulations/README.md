Each line in dataset.details_cases_total.txt and dataset.details_controls_total.txt helps submit a single simulation as a job into nextflow.

dataset.details_cases_total.txt contains below columns.

Column 1 : Integer to identify the founder scenario, values take 1-10   
Column 2 : Integer to track each simulation (1-1320)   
Column 3 : Disease causing variant. We simulated 33 disease causing variants. Disease variants are included in the format of "FTDALS1.chr9.027573483." (DCV.chr.start_bp).    
Column 4 : Expected sharing length in total around the disease variant (0.5,1,2 and 5 in cM)


dataset.details_controls_total.txt contains the same columns as dataset.details_cases_total.txt.   
We only run each simulaion once for controls since controls are not simulated to share with the any other samples. Therefore dataset.details_controls_total.txt contains combinations of columns 1, 2 and 3 only.

We use nextflow to perform the simulation study. Which runs the script, [R/Simulations_FoundHaplo_single_founder_effects.R](https://github.com/bahlolab/FoundHaplo/blob/main/Simulation_scripts/R/Simulations_FoundHaplo_single_founder_effects.R).    
Check [Args_Generate_IBD.R](https://github.com/bahlolab/FoundHaplo/blob/main/Simulation_scripts/R/Args_Generate_IBD.R) and [Simulations_FoundHaplo_single_founder_effects.R](https://github.com/bahlolab/FoundHaplo/blob/main/Simulation_scripts/R/Simulations_FoundHaplo_single_founder_effects.R) to understand parameters.   
   
Specify different path_to_save when running cases and controls.

To run cases,
```bash
module load nextflow
nohup /path/run_simulations/cases/./run_cases.nf 
```

To run controls,
```bash
module load nextflow
nohup /path/run_simulations/cases/./run_controls.nf 
```

Once the FoundHaplo run is completed, each output text file in path_to_save should be concatenated to one text file

Run for cases...
```bash
cat *.txt > /path/results_cases.txt
```bash

Run for controls...
```bash
cat *.txt > /path/results_controls.txt
```
results_cases.txt and results_controls.txt will have columns as below;

1.seed_id  
2.ID  
3.Disease variant  
4.Expected sharing length simulated (0.5,1,2 or 5 cM)  
5.Random left sharing length of the test sample with the simulated disease haplotype (0 for controls)  
6.Random right sharing length of the test sample with the simulated disease haplotype (0 for controls)  
7.Simulated rexp() rate used to splice in founder's genomic region into the disease haplotype in use, this value is an indication of how long the IBD segment of the disease haplotype is.  
8.Random left sharing length of the disease haplotype with the simulated founder   
9.Random right sharing length of the disease haplotype with the simulated founder   
10.Simulated Founder sample ID  
11.Simulated Disease haplotype used to calculate the IBD   
12.Simulated Disease haplotype 1 sample ID  
13.Simulated Disease haplotype 2 sample ID  
14.Simulated Disease haplotype 3 sample ID  
15.Simulated Disease haplotype 4 sample ID  
16.Simulated Disease haplotype 5 sample ID  
17.Disease variant  
18.Simulated Disease haplotype used to calculate the IBD   
19.Test sample ID   
20.Total Log likelihood ratio of IBD  
21.Log likelihood ratio of IBD to the left   
22.Log likelihood ratio of IBD to the right  
23.cM traversed in the HMM from the disease locus (estimated length to the left of the IBD/IBS segment in cM)  
24.cM traversed in the HMM from the disease locus (estimated length to the right of the IBD/IBS segment in cM)  
25.Number of allele mismatches between disease haplotype and the test sample in the HM  
26.Total markers identified in the estimated sharing region  
27.MAF threshold used if any  
28.SNP density in the region per cM  
29.Total number of SNPs in the analysis  
30.Total cM length in the analysis  

31.Type of comparison (simulated_case, truth or control)  
simulated_case : The test sample is being simulated to share the respective disease haplotype.  
truth : The test sample is being simulated to share a different disease haplotype.  
control : The test sample is not simulated to share any disease haplotype i.e its a control sample.

Note: Please remove duplicated columns before analysing the results .txt file.

Now, you can analyse the FoundHaplo simulation results in results_cases.txt and results_controls.txt.
