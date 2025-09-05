#!/usr/bin/bash
# written by Puyam Tondonba Singh

fq=/home/tony/project/reference/Pahang_v4/Musa_acuminata_pahang_v4.fasta

gatk CreateSequenceDictionary -R $fq
