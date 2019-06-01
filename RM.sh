#!/bin/bash
source /net/eichler/vol2/home/mvollger/projects/SDA/env_RM.cfg

fasta=$1
dir=tmpRM


mkdir -p $dir 

RepeatMasker \
        -species human \
        -e wublast \
        -dir $dir \
		-pa $(nproc) \
        $fasta
