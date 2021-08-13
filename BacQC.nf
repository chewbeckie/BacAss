#!/usr/bin/env nextflow

/* 2021-06-09 Xanthomonas Miseq Genome Analysis Workflow
by Johanna Wong */

//-----------------------------------------------------------------------------------------//

/*parameters*/

// runfile with info of sampleId,read1,read2 location
params.runfile// = "$baseDir/test/run_file.csv"

// qcfile with info of sampleId and assembly location
params.qcfile// = "$baseDir/test/qc_file.csv"

params.outputdir// = "$baseDir/test/output"
// blastdb location (need to download from ncbi)
params.blastdb// = "/blastdb/ref_prok_rep_genomes"

//-----------------------------------------------------------------------------------------//

/*Channels*/

// input reads pair
Channel
    .fromPath(params.runfile)
    .splitCsv(header:true)
    .map{ row-> tuple(row.sampleId, file(row.read1), file(row.read2)) }
    .set { samples_ch }

// input assembly
Channel
    .fromPath(params.qcfile)
    .splitCsv(header:true)
    .map{ row-> tuple(row.sampleId, file(row.assemblies)) }
    .into { assembly_ch ; blast_ch ; blobtool_ch}

// join reads and assembly into one channel for mapping
samples_ch.join(assembly_ch).set { bwa_ch }

//-----------------------------------------------------------------------------------------//
/*Processes*/

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



