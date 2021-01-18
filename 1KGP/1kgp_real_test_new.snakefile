import os

reference=config["ref"]
grch38_ref=config.get("grch38_ref", reference)
out_dir=config.get("out_dir", "new_2020")
autosomes=[f"chr{x}" for x in range(1, 23)]
gatk = "java -jar ~/code/aganezov/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar"
# picard = "java -jar ~/code/aganezov/picard-tools-2.4.1/picard.jar"
known_snps_from_dbSNP138 = "/scratch/groups/mschatz1/aganezov/CHM13/variants/1000G/from_ena/resources/Homo_sapiens_assembly38.dbsnp138.vcf"
known_indels_from_mills_1000genomes = "/scratch/groups/mschatz1/aganezov/CHM13/variants/1000G/from_ena/resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
known_indels = "/scratch/groups/mschatz1/aganezov/CHM13/variants/1000G/from_ena/resources/Homo_sapiens_assembly38.known_indels.vcf.gz"


sample = config.get("sample", "ERR3239800")
sample_cram = config.get("cram", os.path.join(out_dir, f"{sample}.final.cram"))
clean_rname = config.get("clean_rnames", True)
enable_gatk_spark = config.get("enable_gatk_spark", False)

rule finish:
    input:
        cram=os.path.join(out_dir, f"{sample}.cram"),
        cram_index=os.path.join(out_dir, f"{sample}.cram.crai"),
        samtools_stats=os.path.join(out_dir, f"{sample}.samtools.stats.txt"),
        mosdepth_stats=os.path.join(out_dir, f"{sample}.mosdepth.global.dist.txt")

rule samtools_stats:
    output: os.path.join(out_dir, "{sample,[^._]+}.samtools.stats.txt"),
    input:
        cram=os.path.join(out_dir, "{sample}.cram"),
        cram_crai=os.path.join(out_dir, "{sample}.cram.crai"),
        target_ref=reference,
    log: os.path.join(out_dir, "log", "{sample}.samtools.stats.txt.log")
    benchmark: os.path.join(out_dir, "benchmarks", "{sample}.samtools_stats.benchmark.txt")
    threads: 24
    shell:
        "samtools stats -r {input.target_ref} --reference {input.target_ref} -@ {threads} {input.cram} > {output} 2> {log}"

rule mosdepth_stats:
    output:
        global_dist=os.path.join(out_dir, "{sample,[^._]+}.mosdepth.global.dist.txt"),
        regions_dist=os.path.join(out_dir, "{sample,[^._]+}.mosdepth.region.dist.txt"),
        summary=os.path.join(out_dir, "{sample,[^._]+}.mosdepth.summary.txt"),
        regions_bed=os.path.join(out_dir, "{sample,[^._]+}.regions.bed.gz"),
        regions_bed_i=os.path.join(out_dir, "{sample,[^._]+}.regions.bed.gz.csi"),
    input:
        cram=os.path.join(out_dir, "{sample}.cram"),
        cram_index=os.path.join(out_dir, "{sample}.cram.crai"),
        chm13_ref=reference,
    log: os.path.join(out_dir, "log", "{sample}.mosdepth.log")
    threads: 24
    benchmark: os.path.join(out_dir, "benchmarks", "{sample}.mosdepth.benchmark.txt")
    params:
        prefix=lambda wc: os.path.join(out_dir, f"{wc.sample}"),
    shell:
        "mosdepth -n --fast-mode --by 500 -t {threads} --fasta {input.chm13_ref} {params.prefix} {input.cram} &> {log}"

rule index_cram:
    output: cram_index=os.path.join(out_dir,"{sample,[^._]+}.cram.crai")
    input: cram=os.path.join(out_dir, "{sample}.cram")
    threads: 24
    shell:
        "samtools index -@ {threads} {input}"

rule bam_to_cram:
    output: os.path.join(out_dir, "{sample,[^._]+}.cram")
    benchmark: os.path.join(out_dir, "benchmark", "{sample}.cram")
    input:
        bam=os.path.join(out_dir, "{sample}.recalibrated.sort.bam"),
        target_ref=reference
    shell: "samtools view -C -T {input.target_ref} -o {output} {input.bam}"

rule recalibrate_bam:
    output: temp(os.path.join(out_dir, "{sample,[^._]+}.recalibrated.sort.bam"))
    input:
        table=os.path.join(out_dir, "{sample}.recal.table"),
        dedup_sorted_bam=os.path.join(out_dir, "{sample}.dedup.sort.bam"),
        dedup_sorted_index=os.path.join(out_dir, "{sample}.dedup.sort.bam.bai"),
        target_ref=reference,
    log: os.path.join(out_dir, "log", "{sample}.recalibrated.sort.bam.log")
    benchmark: os.path.join(out_dir, "benchmark", "{sample}.recalibrated.sort.bam")
    threads: 24
    params:
        gatk_cli="~/code/aganezov/gatk-4.1.9.0/gatk",
        java_options='--java-options \"-Xmx32G -XX:ParallelGCThreads=4 -Djava.io.tmpdir=/dev/shm\"',
        parallelism_options=lambda wc, threads: f"--spark-runner LOCAL --spark-master \'local[{threads}]\'" if enable_gatk_spark else "",
        tool="ApplyBQSRSpark" if enable_gatk_spark else "ApplyBQSR"
    shell:
        "{params.gatk_cli}"
        " {params.tool}"
        " {params.java_options}"
        " {params.parallelism_options}"
        " --preserve-qscores-less-than 6"
        " --static-quantized-quals 10"
        " --static-quantized-quals 20"
        " --static-quantized-quals 30"
        " --tmp-dir /dev/shm"
        " --read-filter GoodCigarReadFilter"
        " -R {input.target_ref}"
        " -O {output}"
        " -I {input.dedup_sorted_bam}"
        " --bqsr-recal-file {input.table}"
        " 2> {log}"

rule get_recalibration_table:
    output: os.path.join(out_dir, "{sample,[^._]+}.recal.table")
    input:
        dedup_sorted_bam=os.path.join(out_dir, "{sample}.dedup.sort.bam"),
        dedup_sorted_index=os.path.join(out_dir, "{sample}.dedup.sort.bam.bai"),
        known_snps_from_dbSNP138=known_snps_from_dbSNP138,
        known_indels=known_indels,
        known_indels_from_mills_1000genomes=known_indels_from_mills_1000genomes,
        target_ref=reference,
    log: os.path.join(out_dir, "log", "{sample}.recal.table.log")
    benchmark: os.path.join(out_dir, "benchmark", "{sample}.recal.table")
    threads: 24
    params:
        autosomes_flag = " ".join(f"-L {c}" for c in autosomes),
        gatk_cli="~/code/aganezov/gatk-4.1.9.0/gatk",
        java_options='--java-options \"-Xmx32G -XX:ParallelGCThreads=4 -Djava.io.tmpdir=/dev/shm\"',
        parallelism_options=lambda wc, threads: f"--spark-runner LOCAL --spark-master \'local[{threads}]\'" if enable_gatk_spark else "",
        tool="BaseRecalibratorSpark" if enable_gatk_spark else "BaseRecalibrator"
    shell:
        "{params.gatk_cli}"
        " {params.tool}"
        " {params.java_options}"
        " {params.parallelism_options}"
        " --preserve-qscores-less-than 6"
        " {params.autosomes_flag}"
        " -R {input.target_ref}"
        " -O {output}"
        " -I {input.dedup_sorted_bam}"
        " --known-sites {input.known_snps_from_dbSNP138}"
        " --known-sites {input.known_indels}"
        " --known-sites {input.known_indels_from_mills_1000genomes}"
        " 2> {log}"

rule index_dedup_bam:
    output: temp(os.path.join(out_dir, "{sample,[^._]+}.dedup.sort.bam.bai"))
    input: os.path.join(out_dir, "{sample}.dedup.sort.bam")
    threads: 24
    log: os.path.join(out_dir, "log", "{sample}.dedup.sort.bam.bai.log")
    shell:
        "samtools index -@ {threads} {input} 2> {log}"



rule mark_duplicates:
    output: temp(os.path.join(out_dir, "{sample,[^._]+}.dedup.sort.bam"))
    input: sorted_bam=os.path.join(out_dir, "{sample}.sort.bam")
    params:
        dedup_metrics = lambda wc: os.path.join(out_dir, f"{wc.sample}.txt"),
        optical_dup_distance=config.get("opt_dup_dis", 100),
        tmp_preifx=lambda wc: os.path.join(out_dir, f"samtools_{wc.sample}", f"{wc.sample}"),
        tmp_dir=lambda wc: os.path.join(out_dir, f"samtools_{wc.sample}"),
    log: os.path.join(out_dir, "log", "{sample}.dedup.sort.bam.log")
    benchmark: os.path.join(out_dir, "benchmark", "{sample}.dedup.bam")
    threads: 24
    shell:
        "mkdir -p {params.tmp_dir} &&"
        " samtools"
        " markdup"
        " -@ {threads}"
        " -S"
        " -T {params.tmp_preifx}"
        " -s"
        " -f {params.dedup_metrics}"
        " -d {params.optical_dup_distance}"
        " {input.sorted_bam}"
        " {output}"
        " 2> {log}"

def aggregate_lanes(wildcards):
    checkpoint_output = os.path.dirname(checkpoints.split_finish.get(**wildcards).output[0])
    return glob_wildcards(os.path.join(
        checkpoint_output,
        f"{wildcards.sample}_1_fastq",
        f"{wildcards.sample}_1." + "{lane}.fastq.gz")
    ).lane

def aggregate_aligned_lanes(wildcards):
    lanes = aggregate_lanes(wildcards)
    return [os.path.join(out_dir, f"{wildcards.sample}.{lane}.sort.bam") for lane in lanes]

rule sorted_merged_bam:
    output: temp(os.path.join(out_dir,"{sample,[^._]+}.sort.bam"))
    input:
        bams=aggregate_aligned_lanes,
        fq1_dir=os.path.join(out_dir, "{sample}_fastq", "{sample}_1_fastq"),
        fq2_dir=os.path.join(out_dir, "{sample}_fastq", "{sample}_2_fastq"),
        f1_split_flag=os.path.join(out_dir, "{sample}_fastq", "done"),
    threads: 24
    benchmark: os.path.join(out_dir,"log","{sample}.nsort.bam.log")
    log: os.path.join(out_dir,"log","{sample}.sort.bam.log")
    shell:
        "samtools merge -@ {threads} {output} {input.bams} 2> {log}"


rule sort_bam_lane:
    output: temp(os.path.join(out_dir,"{sample,[^._]+}.{lane}.sort.bam"))
    input: os.path.join(out_dir,"{sample}.{lane}.fm.bam")
    threads: 5
    benchmark: os.path.join(out_dir,"benchmark","{sample}.{lane}.sort.bam")
    log: os.path.join(out_dir,"log","{sample}.{lane}.sort.bam.log")
    resources:
        time="00:30:00"
    params:
        tmp_prefix=lambda wc: os.path.join(out_dir, f"samtools_tmp_{wc.sample}_{wc.lane}")
    shell:
        "samtools sort -o {output} -@ {threads} -m 1500M {input} -T {params.tmp_prefix} 2> {log}"

rule fixmate_lane:
    output: temp(os.path.join(out_dir,"{sample,[^._]+}.{lane}.fm.bam"))
    input: os.path.join(out_dir,"{sample}.{lane}.bam")
    log: os.path.join("log","{sample,[^_]+}.{lane}.fm.bam.log")
    resources:
        time="2:00:00"
    benchmark: os.path.join(out_dir, "benchmark", "{sample}.{lane}.fm.bam")
    threads: 5
    shell:
        "samtools"
        " fixmate -m"
        " -@ {threads}"
        " {input}"
        " {output}"
        " 2> {log}"

rule align_lane:
    output: temp(os.path.join(out_dir,"{sample,[^._]+}.{lane}.bam"))
    input:
        r1=os.path.join(out_dir, "{sample}_fastq", "{sample}_1_fastq", "{sample}_1.{lane}.fastq.gz"),
        r2=os.path.join(out_dir, "{sample}_fastq", "{sample}_2_fastq", "{sample}_2.{lane}.fastq.gz"),
        target_ref=reference
    threads: 24
    log: os.path.join(out_dir, "log", "{sample}.{lane}.bam.log")
    resources:
        time="3:00:00"
    params:
        rg_tag=lambda
            wc: f"@RG\\tID:{wc.sample}_{wc.lane.replace('.','_')}\\tPL:illumina\\tPM:Unknown\\tLB:{wc.sample}\\tDS:GRCh38\\tSM:{wc.sample}\\tCN:NYGenome\\tPU:{wc.lane.replace('.','_')}",
    benchmark: os.path.join(out_dir, "benchmark", "{sample}.{lane}.bam")
    shell:
        "bwa mem -Y"
        " -K 100000000"
        " -t {threads}"
        " -R \"{params.rg_tag}\""
        " {input.target_ref}"
        " {input.r1}"
        " {input.r2}"
        " 2> {log} | samtools view -Shb -o {output} -"


checkpoint split_finish:
    output: temp(touch(os.path.join(out_dir, "{sample,[^._]+}_fastq", "done")))
    resources:
        time="00:00:05"
    input:
        r1=os.path.join(out_dir, "{sample}_fastq", "{sample}_1_fastq"),
        r2=os.path.join(out_dir,"{sample}_fastq","{sample}_2_fastq"),


rule split_into_lanes:
    output: temp(directory(os.path.join(out_dir, "{sample,[^._]+}_fastq", "{sample}_{n,\d}_fastq"))),
    input: os.path.join(out_dir, "{sample}_{n}.fastq.gz"),
    params:
        sample=lambda wc: wc.sample,
        rn=lambda wc: wc.n,
        sed_command=lambda wc: "sed 's/^@[^ ]* /@/g' |" if clean_rname else "",
        d_threads=lambda wc, threads: max(1, int(threads) // 2),
        c_threads=lambda wc, threads: max(1, int(threads) // 8),
    resources:
        time="3:00:00"
    threads: 24
    benchmark: os.path.join(out_dir, "{sample}_fastq", "log", "{sample}_{n}_fastq")
    shell:
        "mkdir -p {output} && "
        "bgzip -cd -@ {params.d_threads} {input} |"
        " {params.sed_command}"
        " awk "
        "   '"
        "    BEGIN {{FS=\":\"}} "
        "    {{ lane=$3\".\"$4 ; print $0 | \"bgzip -@ {params.c_threads} > {output}/{params.sample}_{params.rn}.\"lane\".fastq.gz\" ;"
        "       for (i = 1; i <= 3; i++) "
        "           {{getline ; print $0 | \"bgzip -@ {params.c_threads} > {output}/{params.sample}_{params.rn}.\"lane\".fastq.gz\"}}"
        "    }}"
        "   '"

rule extract_reads:
    output:
            r1=temp(os.path.join(out_dir, "{sample,[^._]+}_1.fastq.gz")),
            r2=temp(os.path.join(out_dir, "{sample,[^._]+}_2.fastq.gz")),
            z=temp(os.path.join(out_dir, "{sample,[^._]+}_0.fastq.gz")),
            s=temp(os.path.join(out_dir, "{sample,[^._]+}_s.fastq.gz")),
    input:
         source_ref=grch38_ref,
         bam=os.path.join(out_dir, "{sample}.nsorted.bam")
    log: os.path.join(out_dir, "log", "{sample}.reads_extract.log")
    threads: 24
    benchmark: "benchmarks/{sample}.read_extraction.benchmark.txt"
    shell:
        "samtools fastq -@ {threads} -c 2 -1 {output.r1} -2 {output.r2} -0 {output.z} -s {output.s} -n --reference {input.source_ref} {input.bam} &> {log}"

rule name_sort_cram:
    output: temp(os.path.join(out_dir, "{sample_id}.nsorted.bam"))
    input:
        cram=sample_cram,
        source_ref=grch38_ref
    threads: 24
    resources: mem_mb=1500
    params:
        samtools_prefix=lambda wc: os.path.join(out_dir, f"{wc.sample_id}_samtools_nsort_tmp")
    log: os.path.join(out_dir, "log", "{sample_id}.nsorted.cram.log")
    benchmark: os.path.join(out_dir, "benchmarks", "{sample_id}.cram_nsort.benchmark.txt")
    shell:
        "samtools sort -n -@ {threads} -O bam -o {output} -m {resources.mem_mb}M -T {params.samtools_prefix} --reference {input.source_ref} {input.cram} &> {log}"


localrules: split_finish

