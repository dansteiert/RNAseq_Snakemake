rule trim_reads:
    input:
        r1 = os.path.join(work_dir, symlink_dir, "{sample}_R1.fastq.gz"),
        r2 = os.path.join(work_dir, symlink_dir, "{sample}_R2.fastq.gz")
    output:
        r1 = os.path.join(work_dir, trimmed_fastq_dir, "{sample}_R1_trimmed.fastq.gz"),
        r2 = os.path.join(work_dir, trimmed_fastq_dir, "{sample}_R2_trimmed.fastq.gz"),
        html = os.path.join(work_dir, trimmed_fastq_dir, "{sample}_fastp.html"),
        json = os.path.join(work_dir, trimmed_fastq_dir, "{sample}_fastp.json"),
        flag = touch(os.path.join(work_dir, trimmed_fastq_dir, "{sample}.trim.done"))
    params:
        extra = "--trim_poly_x --detect_adapter_for_pe --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 --length_required 20"
    threads: 2
    log:
        os.path.join("logs", "trim_reads", "{sample}.log")
    conda:
        "../envs/fastp.yml"
    shell:
        r"""
        MIN_FASTQ_SIZE=$((1000))  # if you want 1000 MB
            if [ $(stat -Lc%s "{input.r1}") -lt $MIN_FASTQ_SIZE ] || \
               [ $(stat -Lc%s "{input.r2}") -lt $MIN_FASTQ_SIZE ]; then
                echo "⚠️  FASTQ too small for {wildcards.sample}, skipping trimming." | tee -a {log}
                touch {output.r1} {output.r2} {output.html} {output.json}
                exit 0
            fi

        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1} \
            -O {output.r2} \
            --thread {threads} \
            --html {output.html} \
            --json {output.json} \
            {params.extra} \
            2> {log}
        
        # Force filesystem sync for CephFS
        sync
        
        # Wait for files to be fully visible on network filesystem
        sleep 10
        
        # Validate output files are complete and readable
        echo "Validating trimmed outputs..." | tee -a {log}
        
        if [ ! -s {output.r1} ]; then
            echo "ERROR: {output.r1} is empty or missing" | tee -a {log}
            exit 1
        fi
        
        if [ ! -s {output.r2} ]; then
            echo "ERROR: {output.r2} is empty or missing" | tee -a {log}
            exit 1
        fi
        
        # Test gzip integrity
        if ! gunzip -t {output.r1} 2>> {log}; then
            echo "ERROR: {output.r1} is corrupted" | tee -a {log}
            exit 1
        fi
        
        if ! gunzip -t {output.r2} 2>> {log}; then
            echo "ERROR: {output.r2} is corrupted" | tee -a {log}
            exit 1
        fi
        
        echo "✓ Trimming and validation complete for {wildcards.sample}" | tee -a {log}
        """

rule fastqc_trimmed:
    input:
        fastq = os.path.join(work_dir, trimmed_fastq_dir, "{sample}_R{read}_trimmed.fastq.gz"),
        flag = os.path.join(work_dir, trimmed_fastq_dir, "{sample}.trim.done")
    output:
        html = os.path.join(work_dir, "metrics", "trimmed", "{sample}_R{read}_trimmed_fastqc.html"),
        zip = os.path.join(work_dir, "metrics", "trimmed", "{sample}_R{read}_trimmed_fastqc.zip")
    threads: 2
    log:
        os.path.join("logs", "fastqc_trimmed", "{sample}_R{read}.log")
    conda:
        "../envs/fastqc.yml"
    shell:
        """
        mkdir -p {work_dir}/metrics/trimmed
        
        # Additional safety check before processing
        if ! gunzip -t {input.fastq} 2>> {log}; then
            echo "ERROR: Input file {input.fastq} is corrupted" | tee -a {log}
            exit 1
        fi
        
        fastqc -t {threads} -o {work_dir}/metrics/trimmed {input.fastq} 2> {log}
        """