configfile: "config.yaml"

import os
from os.path import basename, join

WD = config['wd']

FASTQ_DIR = join(WD, config['fastq_dir'])
REFGEN_DIR = join(WD, config['refgen_dir'])
DBSNP_DIR = join(WD, config['dbsnp_dir'])

TRIMMED_DIR = join(WD, config['trimmed_dir'])
SAM_DIR = join(WD, config['sam_dir'])
READ_GROUPS_DIR = join(WD, config['read_groups_dir'])
MARK_DUPS_DIR = join(WD, config['mark_dups_dir'])
RECAL_BASES_DIR = join(WD, config['recal_bases_dir'])
VARIANTS_DIR = join(WD, config['variants_dir'])
JOINT_VARIANTS_DIR = join(WD, config['joint_variants_dir'])

WC_fastqs = glob_wildcards(join(FASTQ_DIR, '{sample}.{pair}.fq.gz'))
WC_refgens = glob_wildcards(join(REFGEN_DIR, '{genome}.fa'))

SAMPLES = set(WC_fastqs.sample)
PAIRS = ['R1', 'R2']
GENOMES = set(WC_refgens.genome)


rule all:
    input:
        expand('%s/all.joint.{genome}.vcf' % JOINT_VARIANTS_DIR, genome=GENOMES)
    run:
        print("DNAVCW FINISHED WITH NO EXCEPTIONS!")


rule cut_adapters:
    input:
        r1 = '%s/{sample}.R1.fq.gz' % FASTQ_DIR,
        r2 = '%s/{sample}.R2.fq.gz' % FASTQ_DIR
    output:
        r1 = '%s/{sample}.R1_val_1.fq.gz' % TRIMMED_DIR,
        r2 = '%s/{sample}.R2_val_2.fq.gz' % TRIMMED_DIR
    params:
        outdir = TRIMMED_DIR
    shell:
        'trim_galore -q 30 --illumina --paired {input.r1} {input.r2} -o {params.outdir}'


rule bowtie2_build_genome:
    input:
        "%s/{genome}.fa" % REFGEN_DIR
    output:
        "%s/{genome}.fa.1.bt2" % REFGEN_DIR
    threads: config['bowtie2_build_threads']
    shell:
        "bowtie2-build --threads {threads} -f {input} {input}"


rule bowtie2_map:
    input:
        refgen = '%s/{genome}.fa' % REFGEN_DIR,
        idx = '%s/{genome}.fa.1.bt2' % REFGEN_DIR,
        f1 = '%s/{sample}.R1_val_1.fq.gz' % TRIMMED_DIR,
        f2 = '%s/{sample}.R2_val_2.fq.gz' % TRIMMED_DIR
    output:
        '%s/{sample}.{genome}.sam' % SAM_DIR
    threads: config['bowtie2_map_threads']
    shell:
        "bowtie2 -x {input.refgen} -1 {input.f1} -2 {input.f2} -S {output} -p {threads}"


rule add_read_groups:
    input:
        '%s/{sample}.{genome}.sam' % SAM_DIR
    output:
        '%s/{sample}.{genome}.bam' % READ_GROUPS_DIR
    params:
        picard = config['picard_path']
    shell:
        "{params.picard} AddOrReplaceReadGroups I={input} O={output} SO=coordinate RGID=id "
        "RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=sample; "


rule mark_duplicates:
    input:
        '%s/{sample}.{genome}.bam' % READ_GROUPS_DIR
    output:
        bam = '%s/{sample}.{genome}.bam' % MARK_DUPS_DIR,
        metrics = '%s/{sample}.{genome}.metrics' % MARK_DUPS_DIR
    params:
        picard = config['picard_path']
    shell:
        "{params.picard} MarkDuplicates I={input} O={output.bam} CREATE_INDEX=true "
        "VALIDATION_STRINGENCY=SILENT M={output.metrics}"


rule create_genome_dictionary:
    input:
        "%s/{genome}.fa" % REFGEN_DIR
    output:
        "%s/{genome}.dict" % REFGEN_DIR
    params:
        picard = config['picard_path']
    shell:
        "{params.picard} CreateSequenceDictionary R={input} O={output}"


rule index_genome:
    input:
        "%s/{genome}.fa" % REFGEN_DIR
    output:
        "%s/{genome}.fa.fai" % REFGEN_DIR
    shell:
        "samtools faidx {input}"


rule bgzip_dbsnp_vcf:
    input:
        "%s/{genome}.vcf" % DBSNP_DIR
    output:
        "%s/{genome}.vcf.gz" % DBSNP_DIR
    shell:
        "bgzip {input}"

rule tabix_dbsnp_vcf:
    input:
        "%s/{genome}.vcf.gz" % DBSNP_DIR
    output:
        "%s/{genome}.vcf.gz.tbi" % DBSNP_DIR
    shell:
        "tabix -p vcf {input}"


rule recalibrate_bases:
    input:
        bam = '%s/{sample}.{genome}.bam' % MARK_DUPS_DIR,
        refgen = "%s/{genome}.fa" % REFGEN_DIR,
        refgen_idx = "%s/{genome}.fa.fai" % REFGEN_DIR,
        refgen_dict = "%s/{genome}.dict" % REFGEN_DIR,
        dbsnp = "%s/{genome}.vcf.gz" % DBSNP_DIR,
        dbsnp_idx = "%s/{genome}.vcf.gz.tbi" % DBSNP_DIR
    output:
        '%s/{sample}.{genome}.recal.table' % RECAL_BASES_DIR
    params:
        gatk = config['gatk_path']
    shell:
        "{params.gatk} -T BaseRecalibrator -I {input.bam} -o {output} -R {input.refgen} "
        "-knownSites {input.dbsnp}"

rule print_reads:
    input:
        bam = '%s/{sample}.{genome}.bam' % MARK_DUPS_DIR,
        refgen = "%s/{genome}.fa" % REFGEN_DIR,
        recal = '%s/{sample}.{genome}.recal.table' % RECAL_BASES_DIR
    output:
        '%s/{sample}.{genome}.bam' % RECAL_BASES_DIR
    params:
        gatk = config['gatk_path']
    shell:
        "{params.gatk} -T PrintReads -I {input.bam} -o {output} -R {input.refgen} "
        "-BQSR {input.recal}"


rule variant_calling:
    input:
        bam = '%s/{sample}.{genome}.bam' % RECAL_BASES_DIR,
        refgen = "%s/{genome}.fa" % REFGEN_DIR,
    output:
        "%s/{sample}.{genome}.g.vcf" % VARIANTS_DIR
    params:
        gatk = config['gatk_path']
    shell:
        "{params.gatk} -T HaplotypeCaller -I {input.bam} -o {output} -R {input.refgen} "
        "--genotyping_mode DISCOVERY -stand_call_conf 30 "
        "--emitRefConfidence GVCF"


rule joint_variant_calling:
    input:
        vcfs = expand("%s/{sample}.{genome}.g.vcf" % VARIANTS_DIR, sample=SAMPLES, genome=GENOMES),
        refgen = "%s/{genome}.fa" % REFGEN_DIR
    output:
        '%s/all.joint.{genome}.vcf' % JOINT_VARIANTS_DIR
    run:
        vcfs = '--variant ' + ' --variant '.join(input.vcfs)

        command = "{gatk} -T GenotypeGVCFs -R {refgen} {vcfs} -o {output}".format(gatk=config['gatk_path'],
        refgen=input.refgen, vcfs=vcfs, output=output)

        print(command)
        shell(command)