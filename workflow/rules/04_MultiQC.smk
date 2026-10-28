SAMPLES = config["sample"]

rule multiqc:
    input:
        # Raw FastQC reports
        expand(os.path.join(work_dir, "metrics", "{sample}_R{read}_fastqc.zip"),
               sample=SAMPLES, read=["1", "2"]),
        # Trimmed FastQC reports
        expand(os.path.join(work_dir, "metrics", "trimmed", "{sample}_R{read}_trimmed_fastqc.zip"),
               sample=SAMPLES, read=["1", "2"]),
        # Fastp reports (trimming stats)
        expand(os.path.join(work_dir, trimmed_fastq_dir, "{sample}_fastp.json"),
               sample=SAMPLES)
    output:
        html = os.path.join(work_dir, "multiqc", "multiqc_report.html"),
        data = directory(os.path.join(work_dir, "multiqc", "multiqc_data"))
    params:
        # Search these directories for reports
        search_dirs = f"{work_dir}/metrics {work_dir}/metrics/trimmed {work_dir}/{trimmed_fastq_dir} {work_dir}/{log_dir}",
        outdir = os.path.join(work_dir, "multiqc"),
        title = "RNA-seq QC Report"
    log:
        os.path.join("logs", "multiqc", "multiqc.log")
    conda:
        "../envs/multiqc.yml"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=60,
        nodes=1
    shell:
        """
        mkdir -p {params.outdir}
        multiqc \
            {params.search_dirs} \
            -o {params.outdir} \
            --title "{params.title}" \
            --force \
            --verbose \
            2> {log}
        """