import os

reference = config["ref"]

out_dir = config.get("out_dir","output_var")
alignment_dir = config.get("alignment_dir", "")

gatk_cli = config.get("gatk_cli", "gatk")

known_snps_from_dbSNP138 = config.get("dbsnp38_snps", "/scratch/groups/mschatz1/aganezov/CHM13/variants/1000G/from_ena/resources/Homo_sapiens_assembly38.dbsnp138.vcf")
known_indels_from_mills_1000genomes = config.get("mills_indels", "/scratch/groups/mschatz1/aganezov/CHM13/variants/1000G/from_ena/resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz")
known_indels = config.get("indels", "/scratch/groups/mschatz1/aganezov/CHM13/variants/1000G/from_ena/resources/Homo_sapiens_assembly38.known_indels.vcf.gz")
hapmap = config.get("hapmap", "/scratch/groups/mschatz1/aganezov/CHM13/variants/1000G/from_ena/resources/hapmap_3.3.hg38.vcf.gz")
omni = config.get("omni", "/scratch/groups/mschatz1/aganezov/CHM13/variants/1000G/from_ena/resources/1000G_omni2.5.hg38.vcf.gz")
kg_snps = config.get("kgp_snps", "/scratch/groups/mschatz1/aganezov/CHM13/variants/1000G/from_ena/resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz")

regions = config.get("regions", [f"chr{x}" for x in range(1, 23)])
regions_regex = "(" + "|".join(regions) + ")"

cohort = config["cohort"]
samples = config["samples"]
if not isinstance(samples, list):
    with open(samples, "rt") as source:
        samples = [line.strip() for line in source]

enable_gatk_spark = config.get("enable_gatk_spark", False)

rule finish:
    input: vcf=os.path.join(out_dir,f"{cohort}.recalibrated.snp_indel.vcf")

rule apply_recalibration_indel:
    output: temp(os.path.join(out_dir,"{cohort,[^._]+}.recalibrated.snp_indel.vcf"))
    input:
        vcf=os.path.join(out_dir,"{cohort}.recalibrated.snp.vcf"),
        target_ref=reference,
        recal=os.path.join(out_dir,"{cohort}.recalibrate_indel.recal"),
        tranches=os.path.join(out_dir,"{cohort}.recalibrate_indel.tranches")
    threads: 4
    log: os.path.join(out_dir,"log","{cohort}.apply_vqsr_snp"),
    benchmark: os.path.join(out_dir,"benchmark","{cohort}.apply_vqsr_snp")
    params:
        gatk_cli=gatk_cli,
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
    output: temp(os.path.join(out_dir,"{cohort,[^._]+}.recalibrated.snp.vcf"))
    input:
        vcf=os.path.join(out_dir,"{cohort}.genotyped.vcf"),
        target_ref=reference,
        recal=os.path.join(out_dir,"{cohort}.recalibrate_snp.recal"),
        tranches=os.path.join(out_dir,"{cohort}.recalibrate_snp.tranches"),
    threads: 4
    log: os.path.join(out_dir,"log","{cohort}.apply_vqsr_snp"),
    benchmark: os.path.join(out_dir,"benchmark","{cohort}.apply_vqsr_snp")
    params:
        gatk_cli=gatk_cli,
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
        recal=os.path.join(out_dir,"{cohort,[^._]+}.recalibrate_snp.recal"),
        tranches=os.path.join(out_dir,"{cohort,[^._]+}.recalibrate_snp.tranches"),
        # R=os.path.join(out_dir,"{sample,[^._]+}.recalibrate_snp.plots.R")
    input:
        vcf=os.path.join(out_dir,"{cohort}.genotyped.vcf"),
        target_ref=reference,
        hapmap=hapmap,
        kg_omni=omni,
        kg_snps=kg_snps,
        dbsnp=known_snps_from_dbSNP138,
    threads: 4
    log: os.path.join(out_dir,"log","{cohort}.recalibrate_snp"),
    benchmark: os.path.join(out_dir,"benchmark","{cohort}.recalibrate_snp")
    params:
        gatk_cli=gatk_cli,
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
        recal=os.path.join(out_dir,"{cohort,[^._]+}.recalibrate_indel.recal"),
        tranches=os.path.join(out_dir,"{cohort,[^._]+}.recalibrate_indel.tranches"),
        # R=os.path.join(out_dir,"{sample,[^._]+}.recalibrate_indel.plots.R"),
    input:
        vcf=os.path.join(out_dir,"{cohort}.genotyped.vcf"),
        target_ref=reference,
        dbsnp=known_snps_from_dbSNP138,
        kg_milss=known_indels_from_mills_1000genomes,
    threads: 4
    log: os.path.join(out_dir,"log","{cohort}.recalibrate_indel"),
    benchmark: os.path.join(out_dir,"benchmark","{cohort}.recalibrate_indel")
    params:
        gatk_cli = gatk_cli,
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


rule gather_genotyped_vcfs:
    output: temp(os.path.join(out_dir,"{cohort,[^._]+}.genotyped.vcf"))
    input: [os.path.join(out_dir, "{cohort}." + f"{region}.genotyped.vcf") for region in regions]
    threads: 4
    log: os.path.join(out_dir,"log","{cohort}.genotyped.vcf")
    benchmark: os.path.join(out_dir,"benchmark","{cohort}.genotyped.vcf")
    params:
        gatk_cli = gatk_cli,
        java_options=lambda wc, threads: f'--java-options \"-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads={threads} -Djava.io.tmpdir=/dev/shm\"',
        inputs = lambda wc: " -I ".join(os.path.join(out_dir, f"{wc.cohort}.{region}.genotyped.vcf") for region in regions),
    shell:
        "{params.gatk_cli}"
        " GatherVcfs"
        " -I {params.inputs}"
        " -O {output}"
        " &> {log}"

rule genotype_genomicsdb:
    output: temp(os.path.join(out_dir,"{cohort,[^._]+}.{region}.genotyped.vcf"))
    input:
        target_ref=reference,
        gdb_flag=os.path.join(out_dir, "{cohort}.{region}.gdb.done")
    threads: 4
    log: os.path.join(out_dir,"log","{cohort}.{region}.genotyped.vcf")
    benchmark: os.path.join(out_dir,"benchmark","{cohort}.{region}.genotyped.vcf")
    params:
        gatk_cli = gatk_cli,
        java_options = lambda wc, threads: f'--java-options \"-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads={threads} -Djava.io.tmpdir=/dev/shm\"',
        gdb_dir_path = lambda wc: os.path.join(out_dir, f"{wc.cohort}.{wc.region}.gdb")
    shell:
        "{params.gatk_cli}"
        " GenotypeGVCFs"
        " -R {input.target_ref}"
        " -O {output}"
        " -V gendb://{params.gdb_dir_path}"
        " &> {log}"

rule gvcf2genomicsdb:
    output: os.path.join(out_dir,"{cohort}.{region}.gdb.done")
    input: [os.path.join(out_dir, f"{sample}." + "{region}.hc.vcf") for sample in samples]
    log: os.path.join(out_dir,"log","{cohort}.{region}.gdb")
    benchmark: os.path.join(out_dir,"benchmark","{cohort}.{region}.gdb")
    threads: 4
    params:
        gatk_cli = gatk_cli,
        java_options = lambda wc, threads: f'--java-options \"-Xmx32G -XX:+UseParallelGC -XX:ParallelGCThreads={min(threads,4)} -Djava.io.tmpdir=/dev/shm\"',
        inputs = lambda wc: " -V ".join(os.path.join(out_dir, f"{sample}.{wc.region}.hc.vcf") for sample in samples),
        region = lambda wc: wc.region,
        workspace_flag= lambda wc: "--genomicsdb-workspace-path" if not os.path.exists(os.path.join(out_dir, f"{wc.cohort}.{wc.region}.gdb")) else "--genomicsdb-update-workspace-path",
        output_dir = lambda wc: os.path.join(out_dir, f"{wc.cohort}.{wc.region}.gdb")
    shell:
        "{params.gatk_cli}"
        " {params.java_options}"
        " GenomicsDBImport"
        " -V {params.inputs}"
        " {params.workspace_flag} {params.output_dir}"
        " --tmp-dir /dev/shm"
        " --reader-threads {threads}"
        " -L {params.region}"
        " &> {log}"
        " && touch {output}"

rule haplotype_caller:
    output: temp(os.path.join(out_dir, "{sample,[^._]+}.{chr," + regions_regex + "}.hc.vcf"))
    input:
        target_ref=reference,
        recalibrated_cram=os.path.join(alignment_dir, "{sample}.cram"),
        recalibrated_cram_index=os.path.join(alignment_dir, "{sample}.cram.crai"),
    threads: 1
    log: os.path.join(out_dir,"log","{sample}.{chr}.hc.vcf")
    benchmark: os.path.join(out_dir,"benchmark","{sample}.{chr}.hc.vcf")
    params:
        gatk_cli = gatk_cli,
        java_options = lambda wc, threads: f'--java-options \"-Xmx32G -XX:+UseParallelGC -XX:ParallelGCThreads={min(threads, 4)} -Djava.io.tmpdir=/dev/shm\"',
        parallelism_options = lambda wc, threads: f"--spark-runner LOCAL --spark-master \'local[{threads}]\'" if enable_gatk_spark else "",
        tool = "HaplotypeCallerSpark" if enable_gatk_spark else "HaplotypeCaller",
        region=lambda wc: f"{wc.chr}"
    shell:
        "{params.gatk_cli}"
        " {params.tool}"
        " {params.java_options}"
        " {params.parallelism_options}"
        " -R {input.target_ref}"
        " -I {input.recalibrated_cram}"
        " -L {params.region}"
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



