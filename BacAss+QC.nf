#!/usr/bin/env nextflow

/* 2021-06-09 Xanthomonas Miseq Genome Analysis Workflow
by Johanna Wong */

//-----------------------------------------------------------------------------------------//

/*parameters*/

// runfile
params.runfile = "/home/mbac/JWong/projects_wip/20210609_XanthMiSeqPipeline/test_run.csv"
params.outputdir = "/home/mbac/Desktop/JWong/projects_wip/20210609_XanthMiSeqPipeline/test_output"

// blastdb location
params.blastdb = "/home/mbac/JWong/blastdb/ref_prok_rep_genomes"

//-----------------------------------------------------------------------------------------//

/*Channels*/

// input reads pair
Channel
    .fromPath(params.runfile)
    .splitCsv(header:true)
    .map{ row-> tuple(row.sampleId, file(row.read1), file(row.read2)) }
    .set { samples_ch}

//-----------------------------------------------------------------------------------------//

/*Processes*/

// Step 1 - reads preprocessing
process fastp {
    tag "$sampleId"
    errorStrategy 'ignore'
    maxForks 2
    publishDir "${params.outputdir}/${sampleId}/fastp", mode: 'copy'

    input:
    set sampleId, file(read1), file(read2) from samples_ch

    output:
    set sampleId, 
        file("trimmed_${sampleId}.r1.fq.gz"), 
        file("trimmed_${sampleId}.r2.fq.gz") into trimmed_ch, reads_ch
    file("*.html") into fastp_result

    script:
    """
    fastp -i "${read1}" -I "${read2}"\
      -o "trimmed_${sampleId}.r1.fq.gz"\
      -O "trimmed_${sampleId}.r2.fq.gz"
    mv fastp.html ${sampleId}_fastp.html
    """
}

// Step 2 - assembly
process spades {
    tag "$sampleId"
    errorStrategy 'ignore'
    maxRetries 3
    maxForks 2
    publishDir "${params.outputdir}/${sampleId}/spades", mode: 'copy'
    
    input:
    set sampleId, file(r1), file(r2) from trimmed_ch

    output:
    set sampleId, file("${sampleId}_*.fasta") into quast_ch
    file("*_assembly") into spades_result
    set sampleId, file("${sampleId}_scaffolds.fasta") into assembly_ch

    script:
    """
    spades.py -1 $r1 -2 $r2 \
        -k 21,33,55,77 \
        --careful -t 16 --cov-cutoff auto \
        -o ${sampleId}_assembly
    cp ${sampleId}_assembly/scaffolds.fasta ${sampleId}_scaffolds.fasta
    """
}

//clone assembly channels into multiple channels for QC
assembly_ch.into {
    quast_ch; blast_ch; mapping_ch; blobtool_ch}


// Step 3 - quast assembly quality evaluation
process quast {
    conda 'quast=5.0.2'
    errorStrategy 'ignore'
    tag "$sampleId"
    publishDir "${params.outputdir}/${sampleId}", mode: 'copy'

    input:
    set sampleId, path(assembly) from quast_ch
    
    output:
    file("quast") into quast_result

    script:
    """
    quast.py -o quast $assembly
    """

}

// Step 5 - generate blastn Hits file
process blastn {
    errorStrategy 'ignore'
    tag "$sampleId"
    maxForks 4
    publishDir "${params.outputdir}/$sampleId/blobtools", mode: 'copy'

    input:
    set val(sampleId), file(scaffolds) from blast_ch

    output:
    set val(sampleId), file("${sampleId}.blast.out") into blastn_result
    
    script:

    """
    blastn \
    -task megablast \
    -query $scaffolds \
    -db $params.blastdb \
    -outfmt '6 qseqid staxids bitscore std' \
    -num_threads 8 \
    -max_target_seqs 1 \
    -max_hsps 1 \
    -evalue 1e-25 \
    > ${sampleId}.blast.out
    """

}

// join reads and assembly into one channel for mapping
reads_ch.join(mapping_ch).set { bwa_ch }

// Step 4 - map reads back to assembly
process mapping {
    errorStrategy 'ignore'
    tag "$sampleId"
    maxForks 4
    publishDir "${params.outputdir}/${sampleId}/bamfiles", mode: 'copy'

    input:
    set sampleId, file(r1), file(r2), path(assembly) from bwa_ch

    output:
    set val(sampleId), file("${sampleId}.sorted.bam") into mapping_result

    script:
    """
    bwa index $assembly
    bwa mem -t 8 $assembly $r1 $r2 | \
        samtools view -@8 -bS - | \
        samtools sort -@8 -o ${sampleId}.sorted.bam -
    """

}

// join assembly, blastn result and mapping result together for blobtools
blobtool_ch.join(blastn_result).join(mapping_result).set{ blob_ch }

// Step 5 - blobtools contamination evaluation
process blobtools {
    errorStrategy 'ignore'
    tag "$sampleId"
    publishDir "${params.outputdir}/${sampleId}/blobtools", mode: 'copy'

    input:
    set val(sampleId), 
        file(scaffolds), 
        file(blt),
        file(bam) from blob_ch

    output:
    set file("*.blobDB.json*"), file("*.table.txt"), file("*.bam.cov") into blob_result

    script:
   """
    samtools index $bam
    blobtools create -i $scaffolds -y spades \
        -t $blt -b $bam -x bestsumorder \
        -o ${sampleId}
    blobtools view -i ${sampleId}.blobDB.json \
        -x bestsumorder -r species
    blobtools plot -i ${sampleId}.blobDB.json \
        -x bestsumorder -r genus -l 100
    blobtools plot -i ${sampleId}.blobDB.json \
        -x bestsumorder -r species -l 100
    """

}