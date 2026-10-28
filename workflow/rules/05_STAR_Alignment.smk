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
    resources:
        mem_mb = 64000,
        runtime = 60 * 48,
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
        # FIXED: Use trimmed files instead of original symlinks
        read1 = os.path.join(work_dir, trimmed_fastq_dir, "{sample}_R1_trimmed.fastq.gz"),
        read2 = os.path.join(work_dir, trimmed_fastq_dir, "{sample}_R2_trimmed.fastq.gz"),
        flag = os.path.join(work_dir, trimmed_fastq_dir, "{sample}.trim.done"),  # Added dependency
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
    # FIXED: Use container instead of conda to avoid CPU instruction issues
    container: "docker://quay.io/biocontainers/star:2.7.11b--h43eeafb_0"
    # conda: "../envs/star.yml"  # Commented out - use container instead
    resources:
        mem_mb = 100000,  # INCREASED: More memory for stability
        runtime = 60 * 96,
        nodes = 1,
    threads: 10
    shell:
        """
        # Validate input files before running STAR
        if ! gunzip -t {input.read1} 2>> {log}; then
            echo "ERROR: Input file {input.read1} is corrupted" | tee -a {log}
            exit 1
        fi
        
        if ! gunzip -t {input.read2} 2>> {log}; then
            echo "ERROR: Input file {input.read2} is corrupted" | tee -a {log}
            exit 1
        fi
        
        # Run STAR with additional error handling
        STAR \
            --runThreadN {threads} \
            --runMode alignReads \
            --genomeDir {params.output_dir} \
            --readFilesIn {input.read1} {input.read2} \
            --outFileNamePrefix {params.sample_prefix} \
            --outSAMtype BAM SortedByCoordinate \
            --readFilesCommand zcat \
            --quantMode GeneCounts \
            {params.two_pass} \
            --sjdbGTFfile {input.gtf} \
            --limitBAMsortRAM 50000000000 \
            --outBAMsortingThreadN {threads} \
            &> {log}
        
        # Validate outputs
        if [ ! -s {output.alignment} ]; then
            echo "ERROR: Alignment BAM file is empty or missing" | tee -a {log}
            exit 1
        fi
        
        if [ ! -s {output.read_counts} ]; then
            echo "ERROR: Read counts file is empty or missing" | tee -a {log}
            exit 1
        fi
        
        echo "✓ STAR alignment complete for {wildcards.sample}" | tee -a {log}
        """