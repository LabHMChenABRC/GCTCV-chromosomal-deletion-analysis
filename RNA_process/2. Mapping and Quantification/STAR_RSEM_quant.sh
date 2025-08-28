#!/usr/bin/bash
# written by Hou Bo Han

Fq_In=/home/tony/project/pnas_test
STAR_OUTPUT=/home/tony/project/pnas_test/STAR
RSEM_OUTPUT=/home/tony/project/pnas_test/STAR/RSEM

STAR_REF=/home/public_tools/NAS/share/Reference/Musa_acuminata/STAR_musa_v4_CorrectedMaACA8
RSEM_REF=/home/tony/project/DEG_test/RSEM_ref/corrected/Mav4

[ -d $STAR_OUTPUT ] && rm -r $STAR_OUTPUT
[ -d $RSEM_OUTPUT ] && rm -r $RSEM_OUTPUT
mkdir -p $STAR_OUTPUT $RSEM_OUTPUT

cd $STAR_OUTPUT

for Fq in $Fq_In/*_R1_trimmed.fastq.gz
do
    Filename=$(basename "$Fq")
    Filename=${Filename%%_R1_trimmed.fastq.gz}
    echo "Processing: $Filename"

    STAR --genomeDir $STAR_REF \
         --limitBAMsortRAM 64000000000 \
         --readFilesIn ${Fq_In}/${Filename}_R1_trimmed.fastq.gz ${Fq_In}/${Filename}_R2_trimmed.fastq.gz \
         --readFilesCommand zcat \
         --quantMode TranscriptomeSAM GeneCounts \
         --outSAMtype BAM SortedByCoordinate \
         --outBAMcompression 10 \
         --outWigType bedGraph \
         --outWigStrand Unstranded \
         --outWigNorm RPM \
         --outReadsUnmapped Fastx \
         --alignIntronMax 10000 \
         --outFilterMultimapNmax 1 \
         --outSAMattributes All \
         --twopassMode Basic \
         --outSAMattrRGline ID:$Filename LB:$Filename PL:Illumina PU:$Filename SM:$Filename \
         --runThreadN 24 \
         --outBAMsortingThreadN 16 \
         --outFileNamePrefix "${Filename}_"

    samtools index ${Filename}_Aligned.sortedByCoord.out.bam
    rm -r ${Filename}__STARgenome ${Filename}__STARpass1

    if [ -f ${Filename}_Aligned.toTranscriptome.out.bam ]; then	
        rsem-calculate-expression --paired-end \
                                  --alignments \
                                  --no-bam-output \
                                  --quiet \
                                  -p 24 ${Filename}_Aligned.toTranscriptome.out.bam \
                                  $RSEM_REF \
                                  $RSEM_OUTPUT/$Filename
        rm -r $RSEM_OUTPUT/$Filename.stat
    fi
done

cd $STAR_OUTPUT
STAR_log2table.pl > STAR_Mapping_Info.tsv
