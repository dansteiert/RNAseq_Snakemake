rule fastqc:
    input:
        read1 = os.path.join(work_dir, symlink_dir, "{sample}_R1.fastq.gz"),
        read2 = os.path.join(work_dir, symlink_dir, "{sample}_R2.fastq.gz")
    params:
        outdir = os.path.join(work_dir, "metrics")
    output:
        html = [
            os.path.join(work_dir, "metrics", "{sample}_R1_fastqc.html"),
            os.path.join(work_dir, "metrics", "{sample}_R2_fastqc.html")
        ],
        zip = [
            os.path.join(work_dir, "metrics", "{sample}_R1_fastqc.zip"),
            os.path.join(work_dir, "metrics", "{sample}_R2_fastqc.zip")
        ]
    conda:
        "../envs/fastqc.yml"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=24 * 60,
        nodes=1
    shell:
        "mkdir -p {params.outdir}; fastqc {input.read1} {input.read2} -o {params.outdir}"
