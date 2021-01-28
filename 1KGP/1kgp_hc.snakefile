import os

reference = config["ref"]
grch38_ref = config.get("grch38_ref",reference)
out_dir = config.get("out_dir","new_2020")
autosomes = [f"chr{x}" for x in range(1,23)]
gatk = "java -jar ~/code/aganezov/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar"
known_snps_from_dbSNP138 = "/scratch/groups/mschatz1/aganezov/CHM13/variants/1000G/from_ena/resources/Homo_sapiens_assembly38.dbsnp138.vcf"
known_indels_from_mills_1000genomes = "/scratch/groups/mschatz1/aganezov/CHM13/variants/1000G/from_ena/resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
known_indels = "/scratch/groups/mschatz1/aganezov/CHM13/variants/1000G/from_ena/resources/Homo_sapiens_assembly38.known_indels.vcf.gz"
hapmap = "/scratch/groups/mschatz1/aganezov/CHM13/variants/1000G/from_ena/resources/hapmap_3.3.hg38.vcf.gz"
omni = "/scratch/groups/mschatz1/aganezov/CHM13/variants/1000G/from_ena/resources/1000G_omni2.5.hg38.vcf.gz"
kg_snps = "/scratch/groups/mschatz1/aganezov/CHM13/variants/1000G/from_ena/resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz"

sample = config["sample"]
sample_cram = config.get("cram", os.path.join(out_dir, f"{sample}.final.cram"))
clean_rname = config.get("clean_rnames", True)
enable_gatk_spark = config.get("enable_gatk_spark", False)

rule finish:
    input: vcf=os.path.join(out_dir,f"{sample}.recalibrated.snp_indel.vcf")

rule apply_recalibration_indel:
    output: temp(os.path.join(out_dir,"{sample,[^._]+}.recalibrated.snp_indel.vcf"))
    input:
        vcf=os.path.join(out_dir,"{sample}.recalibrated.snp.vcf"),
        target_ref=reference,
        recal=os.path.join(out_dir,"{sample}.recalibrate_indel.recal"),
        tranches=os.path.join(out_dir,"{sample}.recalibrate_indel.tranches")
    threads: 4
    log: os.path.join(out_dir,"log","{sample}.apply_vqsr_snp"),
    benchmark: os.path.join(out_dir,"benchmark","{sample}.apply_vqsr_snp")
    params:
        gatk_cli="~/code/aganezov/gatk-4.1.9.0/gatk",
        java_options=lambda wc, threads: f'--java-options \"-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads={threads} -Djava.io.tmpdir=/dev/shm\"',
    shell:
        "{params.gatk_cli}"
        " ApplyVQSR"
        " -V {input.vcf}"
        " -O {output}"
        " -mode INDEL"
        " --recal-file {input.recal}"
        " --tranches-file {input.tranches}"
        " -truth-sensitivity-filter-level 99.0"
        " &> {log}"

rule apply_recalibration_snp:
    output: temp(os.path.join(out_dir,"{sample,[^._]+}.recalibrated.snp.vcf"))
    input:
        vcf=os.path.join(out_dir,"{sample}.genotyped.vcf"),
        target_ref=reference,
        recal=os.path.join(out_dir,"{sample}.recalibrate_snp.recal"),
        tranches=os.path.join(out_dir,"{sample}.recalibrate_snp.tranches"),
    threads: 4
    log: os.path.join(out_dir,"log","{sample}.apply_vqsr_snp"),
    benchmark: os.path.join(out_dir,"benchmark","{sample}.apply_vqsr_snp")
    params:
        gatk_cli="~/code/aganezov/gatk-4.1.9.0/gatk",
        java_options= lambda wc, threads: f'--java-options \"-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads={threads} -Djava.io.tmpdir=/dev/shm\"',
    shell:
        "{params.gatk_cli}"
        " ApplyVQSR"
        " -V {input.vcf}"
        " -O {output}"
        " -mode SNP"
        " --recal-file {input.recal}"
        " --tranches-file {input.tranches}"
        " -truth-sensitivity-filter-level 99.8"
        " &> {log}"

rule create_recalibration_snp_data:
    output:
        recal=os.path.join(out_dir,"{sample,[^._]+}.recalibrate_snp.recal"),
        tranches=os.path.join(out_dir,"{sample,[^._]+}.recalibrate_snp.tranches"),
        # R=os.path.join(out_dir,"{sample,[^._]+}.recalibrate_snp.plots.R")
    input:
        vcf=os.path.join(out_dir,"{sample}.genotyped.vcf"),
        target_ref=reference,
        hapmap=hapmap,
        kg_omni=omni,
        kg_snps=kg_snps,
        dbsnp=known_snps_from_dbSNP138,
    threads: 4
    log: os.path.join(out_dir,"log","{sample}.recalibrate_snp"),
    benchmark: os.path.join(out_dir,"benchmark","{sample}.recalibrate_snp")
    params:
        gatk_cli="~/code/aganezov/gatk-4.1.9.0/gatk",
        java_options= lambda wc, threads: f'--java-options \"-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads={threads} -Djava.io.tmpdir=/dev/shm\"',
    shell:
        "{params.gatk_cli}"
        " VariantRecalibrator"
        " {params.java_options}"
        " -R {input.target_ref}"
        " -V {input.vcf}"
        " --mode SNP"
        " -O {output.recal}"
        " --tranches-file {output.tranches}"
        " --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp}"
        " --resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hapmap}"
        " --resource:omni,known=false,training=true,truth=true,prior=12.0 {input.kg_omni}"
        " --resource:1000G,known=false,training=true,truth=false,prior=10.0 {input.kg_snps}"
        " -an QD -an FS -an MQ -an ReadPosRankSum -an MQRankSum -an SOR -an DP"
        " -tranche 100.0 -tranche 99.8 -tranche 99.6 -tranche 99.4 -tranche 99.2 -tranche 99.0 -tranche 95.0 -tranche 90.0"
        " --max-gaussians 6"
        " &> {log}"


rule create_recalibration_indel_data:
    output:
        recal=os.path.join(out_dir,"{sample,[^._]+}.recalibrate_indel.recal"),
        tranches=os.path.join(out_dir,"{sample,[^._]+}.recalibrate_indel.tranches"),
        # R=os.path.join(out_dir,"{sample,[^._]+}.recalibrate_indel.plots.R"),
    input:
        vcf=os.path.join(out_dir,"{sample}.genotyped.vcf"),
        target_ref=reference,
        dbsnp=known_snps_from_dbSNP138,
        kg_milss=known_indels_from_mills_1000genomes,
    threads: 4
    log: os.path.join(out_dir,"log","{sample}.recalibrate_indel"),
    benchmark: os.path.join(out_dir,"benchmark","{sample}.recalibrate_indel")
    params:
        gatk_cli = "~/code/aganezov/gatk-4.1.9.0/gatk",
        java_options = lambda wc, threads: f'--java-options \"-Xmx8G -XX:ParallelGCThreads={threads} -Djava.io.tmpdir=/dev/shm\"',
    shell:
        "{params.gatk_cli}"
        " VariantRecalibrator"
        " {params.java_options}"
        " -R {input.target_ref}"
        " -V {input.vcf}"
        " --mode INDEL"
        " -O {output.recal}"
        " --tranches-file {output.tranches}"
        " --resource:mills,known=false,training=true,truth=true,prior=15.0 {input.kg_milss}"
        " --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp}"
        " -an QD -an FS -an ReadPosRankSum -an MQRankSum -an SOR -an DP"
        " -tranche 100.0 -tranche 99.0 -tranche 95.0 -tranche 92.0 -tranche 90.0"
        " --max-gaussians 4"
        " &> {log}"


rule genotype_vcf:
    output: temp(os.path.join(out_dir,"{sample,[^._]+}.genotyped.vcf"))
    input:
        target_ref=reference,
        gvcf=os.path.join(out_dir, "{sample}.hc.vcf")
    threads: 4
    log: os.path.join(out_dir,"log","{sample}.genotyped.vcf")
    benchmark: os.path.join(out_dir,"benchmark","{sample}.genotyped.vcf")
    params:
        gatk_cli = "~/code/aganezov/gatk-4.1.9.0/gatk",
        java_options = lambda wc, threads: f'--java-options \"-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads={threads} -Djava.io.tmpdir=/dev/shm\"',
    shell:
        "{params.gatk_cli}"
        " GenotypeGVCFs"
        " -R {input.target_ref}"
        " -O {output}"
        " -V {input.gvcf}"
        " &> {log}"

rule haplotype_caller:
    output: os.path.join(out_dir, "{sample,[^._]+}.hc.vcf")
    input:
        target_ref=reference,
        recalibrated_cram=os.path.join(out_dir, "{sample}.cram"),
        recalibrated_cram_index=os.path.join(out_dir, "{sample}.cram.crai"),
    threads: 24
    log: os.path.join(out_dir,"log","{sample}.hc.vcf")
    benchmark: os.path.join(out_dir,"benchmark","{sample}.hc.vcf")
    params:
        gatk_cli = "~/code/aganezov/gatk-4.1.9.0/gatk",
        java_options = lambda wc, threads: f'--java-options \"-Xmx32G -XX:+UseParallelGC -XX:ParallelGCThreads={min(threads, 4)} -Djava.io.tmpdir=/dev/shm\"',
        parallelism_options = lambda wc, threads: f"--spark-runner LOCAL --spark-master \'local[{threads}]\'" if enable_gatk_spark else "",
        tool = "HaplotypeCallerSpark" if enable_gatk_spark else "HaplotypeCaller"
    shell:
        "{params.gatk_cli}"
        " {params.tool}"
        " {params.java_options}"
        " {params.parallelism_options}"
        " -R {input.target_ref}"
        " -I {input.recalibrated_cram}"
        " -pairHMM AVX_LOGLESS_CACHING"
        " -O {output}"
        " -ERC GVCF"
        " -A DepthPerAlleleBySample"
        " -A DepthPerSampleHC"
        " -A InbreedingCoeff"
        " -A StrandBiasBySample"
        " -A Coverage"
        " -A MappingQualityRankSumTest"
        " -A MappingQualityZero"
        " -A QualByDepth"
        " -A RMSMappingQuality"
        " -A ReadPosRankSumTest"
        " &> {log}"



