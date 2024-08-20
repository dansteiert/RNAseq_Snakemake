rule create_index:
    input:
        genome = genome,
        gtf = gtf
    output:
        genome = os.path.join(genome_index_dir, f"RNA_STAR_{genome_version}_length_{read_length}", "Genome"),
        SA = os.path.join(genome_index_dir, f"RNA_STAR_{genome_version}_length_{read_length}", "SA"),
        SAindex = os.path.join(genome_index_dir, f"RNA_STAR_{genome_version}_length_{read_length}", "SAindex"),
    params:
        read_length = read_length,
        output_dir = os.path.join(genome_index_dir, f"RNA_STAR_{genome_version}_length_{read_length}"),
        star_version = star_version,
    log:
        os.path.join(work_dir, log_dir, "create_index.log")
    # container: "docker://mgibio/star"
    conda: "../envs/star.yml"
    resources:
        mem_mb = 64000,
        runtime = 60 * 24,
        nodes=1,
    threads: 10
    shell:
        "STAR "
        "--runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir {params.output_dir} "
        "--genomeFastaFiles {input.genome} "
        "--sjdbGTFfile  {input.gtf} "
        "--sjdbOverhang {params.read_length} "
        "&> {log}"


rule star_alignment:
    input:
        genome = os.path.join(genome_index_dir, f"RNA_STAR_{genome_version}_length_{read_length}", "Genome"),
        SA = os.path.join(genome_index_dir, f"RNA_STAR_{genome_version}_length_{read_length}", "SA"),
        SAindex = os.path.join(genome_index_dir, f"RNA_STAR_{genome_version}_length_{read_length}", "SAindex"),
        read1 = os.path.join(work_dir, symlink_dir, "{sample}_R1.fastq.gz"),
        read2 = os.path.join(work_dir, symlink_dir, "{sample}_R2.fastq.gz"),
        gtf = gtf,
    output:
        alignment = os.path.join(work_dir, alignment_dir, "{sample}_Aligned.sortedByCoord.out.bam"),
        read_counts = os.path.join(work_dir, alignment_dir, "{sample}_ReadsPerGene.out.tab"),
    params:
        read_length = read_length,
        output_dir = os.path.join(genome_index_dir, f"RNA_STAR_{genome_version}_length_{read_length}"),
        star_version = star_version,
        sample_prefix = lambda wildcards: os.path.join(work_dir, alignment_dir, f"{wildcards.sample}_"),
        two_pass = "--twopassMode Basic " if two_pass else "",
    log:
        os.path.join(work_dir, log_dir, "alignment_{sample}.log")
    # container: "docker://mgibio/star"
    conda: "../envs/star.yml"
    resources:
        mem_mb = 64000,
        runtime = 60 * 24,
        nodes = 1,
    threads: 10
    shell:
        "STAR "
        "--runThreadN {threads} "
        "--runMode alignReads "
        "--genomeDir {params.output_dir} "
        "--readFilesIn {input.read1} {input.read2} "
        "--outFileNamePrefix {params.sample_prefix} "
        "--outSAMtype BAM SortedByCoordinate "
        "--readFilesCommand zcat "
        "--quantMode GeneCounts "
        "{params.two_pass} "
        "--sjdbGTFfile {input.gtf} "
        "&> {log}"