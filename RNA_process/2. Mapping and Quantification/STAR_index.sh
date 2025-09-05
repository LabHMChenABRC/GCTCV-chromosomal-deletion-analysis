#!/usr/bin/bash
# written by Hou Bo Han

Reference=/home/public_tools/NAS/share/Reference/Musa_acuminata/Musa_acuminata_DH-Pahang_v4.3/Musa_acuminata_pahang_v4.genome.fasta
GFF=/home/public_tools/NAS/share/Reference/Musa_acuminata/Musa_acuminata_DH-Pahang_v4.3/Musa_acuminata_pahang_v4.CorrectedMaACA8.gtf
genomeDir=/home/public_tools/NAS/share/Reference/Musa_acuminata/STAR_musa_v4_CorrectedMaACA8
# MusaV4.3 contains incorrected gene model Macma4_05_g08140.1
# It's replaced with MusaV2 gene models Ma05_g07830(Ma05_t07830.1) and Ma05_g07840(Ma05_t07840.1) 
[ -d $genomeDir ] && rm -r $genomeDir
mkdir -p $genomeDir
STAR \
--runThreadN 32 \
--runMode genomeGenerate \
--genomeDir $genomeDir \
--genomeFastaFiles $Reference \
--genomeSAindexNbases 13 \
--sjdbGTFfile $GFF \
--sjdbGTFtagExonParentTranscript transcript_id \
--sjdbGTFtagExonParentGene gene_id
