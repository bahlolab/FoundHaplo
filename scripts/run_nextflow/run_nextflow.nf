#!/usr/bin/env nextflow

Channel
.fromPath('PATH/manifest.txt')
.splitText(by: 1)
.set{ chunks_ch }

process package_test {

memory { 5.GB * task.attempt }
errorStrategy { task.exitStatus == 143 ? 'ignore' : 'retry' }
maxRetries 4
cpus 1
time '5h'


module 'R/4.2.0'

input:
val x from chunks_ch

output:
stdout result


"""


Rscript /PATH/FoundHaplo/scripts/Args_Generate_FH_score.R $x


"""
}

result.view { it.trim() }
