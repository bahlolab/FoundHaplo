#!/usr/bin/env nextflow

Channel
.fromPath('/mypath/manifest.txt') // edit : give the correct the path
.splitText(by: 1)
.set{ chunks_ch }

process Run_FoundHaplo {

memory { 5.GB * task.attempt } // edit : change the memory requirement as necessary
errorStrategy { task.exitStatus == 143 ? 'ignore' : 'retry' }
maxRetries 4
cpus 1
time '5h' // edit : Change the time requirement as necessary


module 'R/4.2.0' // edit : load the R version with FoundHaplo installed

input:
val x from chunks_ch

output:
stdout result


"""

// edit : give the full path to FoundHaplo/scripts/run_nextflow/Args_Generate_FH_score.R
Rscript /mypath/FoundHaplo/scripts/run_nextflow/Args_Generate_FH_score.R $x


"""
}

result.view { it.trim() }
