#!/usr/bin/env nextflow

// Test Github push 201029
// redo test github blablabla
//params.reads = "/data/mmaj/data/pipelinetests/nextflow2/*.R{1,2}.fq.gz"


// cmd line parameters:
params.data = null
params.genome = "/data/mmaj/genomes/b37/human_g1k_v37.fasta"
params.outdir = 'results'
params.bwaindex = "/data/mmaj/genomes/b37/human_g1k_v37.fasta"
params.rundir= "${launchDir.baseName}" // get basename of dir where script is started
params.bed_WGS = "/data/mmaj/genomes/b37/interval.files/IDT.genetikpanel.GV3/IDT.GV3.ROI.bed"
params.dbsnp= "/data/mmaj/genomes/share.hg19.b37/dbsnp147/All_20160601.vcf"
params.gatkTEMP="${launchDir.baseName}/gatkTEMP"
params.qualimap_bed = null

outdir_full_path= "${launchDir}/${params.outdir}/"
// If datafolder is selected using --data at cmdline, set params.reads to new filelocation:
if (params.data) {
    log.info "manually set input data folder to $params.data"
    params.reads="${params.data}/*R{1,2}*{fq,fastq}.gz"
}
else if (!params.data) {
    log.info "input folder not set, using data in current folder."
    params.reads = "${launchDir}/*R{1,2}*{fq,fastq}.gz"
}

ch_qualimap_bed = params.qualimap_bed ? Channel.value(file(params.qualimap_bed)) : ''


// program related paths (executables, configs etc.):
snpeff="/data/mmaj/programmer/snpeff43/snpEff/snpEff.jar"
manta_config="/data/mmaj/programmer/manta-1.6.0/bin/configManta.py"
multiqc_config="/data/mmaj/programmer/MultiQC/multiqc_config_MJ.yaml"
vcfanno_path="/data/mmaj/genomes/b37/vcfanno.stable.annotations"
snpeff_config="/data/mmaj/programmer/snpeff43/snpEff/snpEff.config"
expansionhunter_b37_db="/data/mmaj/programmer/ExpansionHunter.v4.0.1/variant_catalog/grch37/variant_catalog.json"


// path to file variables:
KGindels="/data/mmaj/genomes/b37/databases/1000G_phase1.indels.b37.vcf"
KGmills="/data/mmaj/genomes/b37/databases/Mills_and_1000G_gold_standard.indels.b37.vcf"
KG_p1_High_snps="/data/mmaj/genomes/b37/databases/1000G_phase1.snps.high_confidence.b37.vcf"
SEref="GRCh37.75"
hapmap="/data/mmaj/genomes/b37/databases/hapmap_3.3.b37.vcf"
omni="/data/mmaj/genomes/b37/databases/1000G_omni2.5.b37.vcf"

EV7_roi = "/data/mmaj/genomes/b37/interval.files/exomes/IDT.exomes.EV7/EV7.ROI.bed"
GV3_roi = "/data/mmaj/genomes/b37/interval.files/IDT.genetikpanel.GV3/IDT.GV3.ROI.bed"
CV4_roi = "/data/mmaj/genomes/b37//interval.files/IDT/IDT.cancerpanel.v1.ROI.bed"

// target region files:
excluderegions="/data/mmaj/genomes/b37/interval.files/whole_genome/gaplist.b37.cleaned.bed"


//GATK_TARGET="-XL ${excluderegions}"

GATK_TARGET="-L ${CV4_roi}"

def summary = [:]
summary['Input']             = params.reads
summary['Qualimap_bed']      = params.qualimap_bed
summary['Genome']            = params.genome



runID="TESTRUNID"
user="$USER"
runtype="testNF"


log.info """\
=======================================
KGA Vejle WGS nextflow pipeline v1
=======================================
genome   : $params.genome
datafolder: $params.data
reads    : $params.reads
results  : $params.outdir
user     : $user
rundir   : $params.rundir
GATK target: $GATK_TARGET

"""

Channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs_ch }

process trim_galore_std {
    cpus 20
    tag "$sampleID"

    input:
    tuple val(sampleID), path(reads) from read_pairs_ch

    output:
    tuple sampleID, path("*_val_*") into trimmed_reads

    script:
    """
    trim_galore $reads --paired
    """
}

process bwa_mem_samblaster{
    tag "$sampleID"
    cpus 10

    input:
    tuple sampleID, path(reads) from trimmed_reads
    output:
    tuple sampleID, path("${sampleID}.samblaster.bam") into samblaster_bam_ch
    shell:
    '''
    bwa0715 mem !{params.bwaindex} -R "@RG\tID:!{runID}\tLB:!{runtype}\tPL:illumina\tSM:!{sampleID}" -t !{task.cpus} !{reads} | samblaster | sambamba view -t 10 -S -f bam /dev/stdin | sambamba sort -t 10 --tmpdir=/data/TMP.!{user}/ -o !{sampleID}.samblaster.bam /dev/stdin
    sambamba index !{sampleID}.samblaster.bam

    '''
}


process baserecalGATK_spark {
    cpus 10
    publishDir "${params.outdir}/${sampleID}/BAM/", mode: 'copy'
    tag "$sampleID"
    
    input:
    tuple sampleID, path(bam) from samblaster_bam_ch
    output:
    tuple sampleID, path("${sampleID}.final.bam") into (varcall_input,bamtools_input, qualimap_input, samtools_input, fastqc_input, cnvnator_input, expansionhunter_input)
    tuple sampleID, path("${sampleID}.*.bai") into expansionhunter_index
   
    shell:
    '''
    gatk419 BQSRPipelineSpark --spark-master local[!{task.cpus}] -R !{params.genome} -I !{bam} --verbosity WARNING --spark-verbosity WARN --known-sites !{KGmills} --known-sites !{KGindels} -O !{sampleID}.final.bam
    sambamba index !{sampleID}.final.bam
    '''
}

process GATK_haplotypecaller{
    cpus 10
    tag "$sampleID"
    publishDir "${params.outdir}/${sampleID}/variantcalls/", mode: 'copy', pattern: "*.HC.*"
    publishDir "${params.outdir}/gvcf/", mode: 'copy', pattern: "*.g.*"
    input:
    tuple sampleID, path(bam) from varcall_input
    
    output:
    tuple sampleID, path("${sampleID}.g.vcf") into sample_gvcf_tuple
    path("${sampleID}.g.vcf") into sample_gvcf_list
    //val true into HC_done_ch
    file "${sampleID}.HC.*"
    file "${sampleID}.g.*"
    file "${sampleID}.HCbamout.*"

    shell:
    '''
     gatk419 HaplotypeCaller -I !{bam} -O !{sampleID}.g.vcf -R !{params.genome} -ERC GVCF !{GATK_TARGET} -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation --dont-use-soft-clipped-bases -bamout !{sampleID}.HCbamout.bam
     gatk419 GenotypeGVCFs -R !{params.genome} -V !{sampleID}.g.vcf -O !{sampleID}.HC.vcf -G StandardAnnotation -G AS_StandardAnnotation
    '''
}

sample_gvcf_list
    .map{" -V "+it}
    .set{gvcflist_done}

gvcflist_done
    .collectFile(name: "${params.outdir}/collectfileTEST1.txt", newLine: false)
    .map {it.text.trim()}.set {gvcfsamples_for_GATK}

process GATK_jointgenotyping{
    cpus 10
    publishDir "${params.outdir}/merged_data/", mode: 'copy'
    input:
    val x from gvcfsamples_for_GATK
    //tag true from HC_done_ch
    output:
    //tuple params.rundir, path("${params.rundir}.merged.g.vcf") into merged_vcf
    path("${params.rundir}.merged.g.vcf") into merged_gVCF
    path("${params.rundir}.merged.RAW.vcf") into merged_RAW_vcf
    shell:
    '''
    gatk419 CombineGVCFs -R !{params.genome} !{x} -O !{params.rundir}.merged.g.vcf !{GATK_TARGET} -G StandardAnnotation -G AS_StandardAnnotation 

    gatk419 GenotypeGVCFs -R !{params.genome} -V !{params.rundir}.merged.g.vcf \
    -O !{params.rundir}.merged.RAW.vcf  \
    !{GATK_TARGET} \
    -G StandardAnnotation -G AS_StandardAnnotation -A SampleList \
    -D !{params.dbsnp}
    '''     
}

process bamtools_STATS {
    publishDir "${params.outdir}/${sampleID}/QC/", mode: 'copy'

    input:
    tuple sampleID, path(bam) from bamtools_input
    output:
    path("${sampleID}.bamtools.sample.stats.txt") into bamtools_ch

    shell:
    '''
    bamtools stats -in !{bam} -insert > !{sampleID}.bamtools.sample.stats.txt
    '''
}

process samtools_STATS {
    publishDir "${params.outdir}/${sampleID}/QC/", mode: 'copy'

    input:
    tuple sampleID, path(bam) from samtools_input
    output:
    path("${sampleID}.samtools.sample.stats.txt") into samtools_ch

    shell:
    '''
    samtools stats !{bam} > !{sampleID}.samtools.sample.stats.txt
    '''
}

process qualimap_STATS {
    publishDir "${params.outdir}/${sampleID}/QC/bamQC", mode: 'copy'

    input:
    tuple sampleID, path(bam) from qualimap_input
    path(targetBED) from ch_qualimap_bed
    output:
    path ("${bam.baseName}/") into bamQCReport
    //path ("*_results.txt") into bamQCReport

    script:
    use_bed = params.qualimap_bed ? "-gff ${targetBED}" : ''
    """
    qualimap bamqc -outdir ${bam.baseName} -bam ${bam} ${use_bed}
    """
}
//    qualimap bamqc -outdir !{bam.baseName} -bam !{bam} -gff !{targetBED}


bamQCReport = bamQCReport.dump(tag:'BamQC')

process fastqc_bam {
    cpus 1
    publishDir "${params.outdir}/${sampleID}/QC/", mode: 'copy',saveAs: { filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename" }
    input:
    tuple sampleID, path(bam) from fastqc_input
    
    output:
    path "*_fastqc.{zip,html}" into fastqc_BAM_results
    shell:
    '''
    fastqc --quiet --threads !{task.cpus} !{bam}
    '''
}

process multiQC {
    publishDir "${params.outdir}/", mode: 'copy'

    input:
    path("*_fastqc.*") from fastqc_BAM_results.collect()
    path("${sampleID}.samtools.sample.stats.txt") from samtools_ch.collect()
    path("${sampleID}.bamtools.sample.stats.txt") from bamtools_ch.collect()
    path("bamQC/*") from bamQCReport.collect()
    output:
    path ("*multiqc_report.html")
    script:
    """
    multiqc -f -q -c ${multiqc_config} --title TESTOMG22 ${outdir_full_path}/*/QC/
    """

}

// Structural variant calling: Manta, CNVnator, ExpansionHunter
/*
process CNVnator{
    publishDir "${params.outdir}/${sampleID}/StructuralVariants/CNVnator", mode: 'copy'
    input:
    tuple sampleID, path(bam) from cnvnator_input

    output:
    path "${sampleID}.CNVnator.output.txt" into cnvnator_out_ch
    shell:
    '''
    source /data/mmaj/programmer/ROOT/bin/thisroot.sh
    cnvnator -root !{sampleID}.CNVnator.root -genome !{params.genome} -tree !{bam}
    cnvnator -root !{sampleID}.CNVnator.root -his 1000 -fasta !{params.genome} 
    cnvnator -root !{sampleID}.CNVnator.root -stat 1000
    cnvnator -root !{sampleID}.CNVnator.root -partition 1000
    cnvnator -root !{sampleID}.CNVnator.root -call 1000 > !{sampleID}.CNVnator.output.txt
    '''
}
*/

process ExpansionHunter{
    publishDir "${params.outdir}/${sampleID}/StructuralVariants/expansionHunter", mode: 'copy'
    input:
    tuple sampleID, path(bam), path(bai) from expansionhunter_input.join(expansionhunter_index)

    output:
    path "${sampleID}.expansionhunter*" into expansionhunter_out_ch

    shell:
    '''
    ExpansionHunter --reads !{bam} --reference !{params.genome} --variant-catalog !{expansionhunter_b37_db} --output-prefix !{sampleID}.expansionhunter
    '''
}









/*
    gatk419 --java-options "-Dsamjdk.buffer_size=131072 -Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx16g" BQSRPipelineSpark --spark-master local[12] -R !{params.genome} -I !{bam} --verbosity WARNING --spark-verbosity WARN --known-sites !{KGmills} --known-sites !{KGindels} -O !{sampleID}.MD.BQSR_spark.bam
*/














/*

sambamba chunk:

    gatk412 BaseRecalibrator -I b4.bam -R !{params.genome} --known-sites !{KGmills} --known-sites !{KGindels} -O gatk1.intervals -L !{params.bed}
    gatk412 ApplyBQSR -R !{params.genome} -I b4.bam --bqsr-recal-file gatk1.intervals -O !{sampleID}.bam
    rm b[1-6]*
    sambamba index !{sampleID}.bam




    gatk412 BaseRecalibrator -I b4.bam -R !{params.genome} --known-sites !{KGmills} --known-sites !{KGindels} -O gatk1.intervals -L !{params.bed}
    gatk412 ApplyBQSR -R !{params.genome} -I b4.bam --bqsr-recal-file gatk1.intervals -O !{sampleID}.bam
    rm b[1-6]*
    sambamba index !{sampleID}.bam


### DEPRECATED COMPLETE PROCESSSES:

process baserecalGATK {
    cpus 20
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${sampleID}/BAM/", mode: 'copy'
    tag "$sampleID"
    
    input:
    tuple sampleID, path(bam) from intermediatebam1
    output:
    tuple sampleID, path("${sampleID}.MD.BQSR.bam") into finalbam_BQSR
    shell:
    '''
    gatk419 BaseRecalibrator -R !{params.genome} -I !{bam} --known-sites !{KGmills} --known-sites !{KGindels} -O gatk1.intervals !{GATK_TARGET}
    gatk419 ApplyBQSR -R !{params.genome} -I !{bam} --bqsr-recal-file gatk1.intervals -O !{sampleID}.MD.BQSR.bam
    sambamba index !{sampleID}.MD.BQSR.bam
    '''
}

process bwa_mem_sambamba{
    //publishDir "${params.outdir}/${sampleID}/BAM/sambamba/", mode: 'copy'
    tag "$sampleID"
    cpus 50

    input:
    tuple sampleID, path(reads) from trimmed_reads
    //tuple sampleID, path(reads) from read_pairs_ch

    
    output:
    tuple sampleID, path("*.bam") into finalbam_sambamba
    //path "*.bai" 
    shell:
    '''
    bwa0715 mem !{params.bwaindex} -R "@RG\tID:!{runID}\tLB:!{runtype}\tPL:illumina\tSM:!{sampleID}" -t !{task.cpus} !{reads} | sambamba view -t 10 -S -f bam /dev/stdin | sambamba sort -t 10 --tmpdir=/data/TMP.!{user}/ -o b2.bam /dev/stdin 

    sambamba markdup --tmpdir=/data/TMP.!{user}/ -t 2 b2.bam b4.bam
    sambamba index b4.bam
    rm b2.*
    '''
}










*/
