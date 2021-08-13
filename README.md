# BacAss
a generic bacteria genome assembly and qc nextflow pipeline

The nextflow pipeline will use illumina short reads (fastq files) for bacterial isolate genomic DNA as input and run the following processes:
1. reads qc by [fastp](https://github.com/OpenGene/fastp)
2. genome assembly by [spades](https://github.com/ablab/spades)
3. genome qc by [quast](https://github.com/ablab/quast)
4. read mapping back to the assembled genome using [bwa](https://github.com/lh3/bwa)
5. blastn check on scaffold basing on the refseq genome database for prokaryotes
6. finally, a qc and taxonomic partitioning on the genome with [blobtools](https://github.com/DRL/blobtools).

## Preparation

### install required tools
This pipeline is dependent on **conda** and **nextflow**. A new environment on conda need to be created by following command:
```
conda env create -f BacAss.yml
```
Once that's done, edit the `nextflow.config` file to specify the location of `process.conda` to the location of `BacAss` conda environment.

### download blast database
The [BLAST database](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) is required for blastn to be ran locally. The blobtools recommend using the `nt` database but it is rather large, so I am using the `ref_prok_rep_genomes` which contains >5700 prokaryotic refseq representative genomes instead. One can download the database by running the following script that come with `blast+` package:
```
update_blastdp.pl --decompress ref_prok_rep_genomes
```

### download NCBI taxdump and create nodesdb for blobtools
This is explained in the [blobtools](https://github.com/DRL/blobtools) github site. Run the following:
```
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -P data/
tar zxf data/taxdump.tar.gz -C data/ nodes.dmp names.dmp
blobtools nodesdb --nodes data/nodes.dmp --names data/names.dmp
```

## Usage
1. Prepare a `csv` run-file that contains the information of the sample name, path to read 1 and read 2. The header must be set as `sampleId`, `read1`, `read2`. Here's an example:
```
sampleId,read1,read2
s1,path/to/s1-read1.fastq.gz,path/to/s1-read2.fastq.gz
s2,path/to/s2-read1.fastq.gz,path/to/s2-read2.fastq.gz
```

2. Run BacAss.nf and define output location and blastdb location, I like to also activate the report function to report the process
```
./BacAss.nf --runfile runfile.csv --outputdir output_path -with-report log.html
```

3. After the above completed, generate a `csv` qc-file that contains the information of sampleId and assembly location
```
sampleId,assemblies
s1,path/to/generated/assembly.fasta
```

4. Run BacQC.nf with the qc-file, run-file and blastdb location specified:
```
./BacQC.nf --runfile runfile.csv \
--qcfile qcfile.csv \
--outputdir outputdir_path \
--blastdb path/to/ref_prok_rep_genomes \
with-report log.html
```