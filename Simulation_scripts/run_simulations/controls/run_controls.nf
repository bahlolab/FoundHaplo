#!/usr/bin/env nextflow

gen_error=1
S_Er=1
IBD_version=1.1
gen_error.rate=0.01
S.Er.rate=20.05
path_to_save="/path/"
root="/path/"
path_MIS_ref="/path/"
path_hapmap="/path/"
path_gnomad_frq="/path/"



meiosis=1

Channel
.fromPath('/path/dataset.details_controls_total.txt')
.splitText(by: 1)
.set{ chunks_ch }




process reverse {


memory { 20.GB * task.attempt }
errorStrategy { task.exitStatus == 143 ? 'ignore' : 'retry' }
maxRetries 4
cpus 1
time '48h'


module 'R/3.6.1'

input:
val x from chunks_ch

output:
stdout result


"""

export TMPDIR=/path/tmp

Rscript /path/Args_Generate_IBD.R $gen_error $S_Er $IBD_version $gen_error.rate $S.Er.rate $path_to_save $meiosis $root $path_MIS_ref $path_hapmap $path_gnomad_frq $x



"""
}

result.view { it.trim() }

