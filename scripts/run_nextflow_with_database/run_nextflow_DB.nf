#!/usr/bin/env nextflow

Channel
.fromPath('/mypath/manifest.txt') // edit
.splitText(by: 1)
.set{ chunks_ch }

process Run_FoundHaplo {

memory { 5.GB * task.attempt } // edit
errorStrategy { task.exitStatus == 143 ? 'ignore' : 'retry' }
maxRetries 4
cpus 1
time '5h' // edit


module 'R/4.2.0' // edit

input:
val x from chunks_ch

output:
stdout result


"""

echo "edit below line with the correct path"
Rscript /mypath/FoundHaplo/scripts/run_nextflow/Args_Generate_FH_score.R $x


"""
}

result.view { it.trim() }
