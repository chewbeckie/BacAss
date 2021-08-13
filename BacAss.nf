#!/usr/bin/env nextflow

/* 2021-06-09 Xanthomonas Miseq Genome Analysis Workflow
by Johanna Wong */

//-----------------------------------------------------------------------------------------//

/*parameters*/

// runfile
params.runfile// = "$baseDir/test/run_file.csv"
// output location
params.outputdir// = "$baseDir/test/output.dir"

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
    cp ${sampleId}_assembly/contigs.fasta ${sampleId}_contigs.fasta
    cp ${sampleId}_assembly/*/final_contigs.fasta ${sampleId}_finalcontigs.fasta
    """
}

//clone assembly channels into multiple channels for QC
assembly_ch.into {
    quast_ch; mapping_ch; blast_ch; blob_ch
}


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

mapping_ch.join(reads_ch).set {mapping_reads_ch}

// Step 4 - map reads back to assembly
process mapping {
    errorStrategy 'ignore'
    tag "$sampleId"
    publishDir "${params.outputdir}/${sampleId}/bamfiles", mode: 'copy'

    input:
    set sampleId, path(assembly), file(r1), file(r2) from mapping_reads_ch

    output:
    file("*.bam") into mapping_result

    script:
    """
    bwa index $assembly
    bwa mem $assembly $r1 $r2 | \
        samtools view -bS - | \
        samtools sort -o ${sampleId}.sorted.bam -
    """

}

/*// Step 5 - generate blastn Hits file
process blastn {
    tag "$sampleId"
    
    input:
    set sampleId, file(scaffolds) from blast_ch

    output:
    file("*.blast.out") into blastn_result
    
 
    script:
    """
    blastn \
    -task megablast \
    -query $scaffolds \
    -db nt -remote \
    -outfmt '6 qseqid staxids bitscore std' \
    -max_target_seqs 10 \
    -max_hsps 1 \
    -evalue 1e-25 \
    > ${scaffolds}.blast.out
    """

}

// Step 5 - blobtools contamination evaluation
process blobtools {
    errorStrategy 'ignore'
    tag "$sampleId"
    publishDir "${params.outputdir}/${sampleId}/blobtools", mode: 'copy'

    input:
    set sampleId, file(scaffolds) from blob_ch
    file(blt) from blastn_result
    file(bam) from mapping_result

    output:
    file("*.blobDB.json*") into blob_result

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
*/