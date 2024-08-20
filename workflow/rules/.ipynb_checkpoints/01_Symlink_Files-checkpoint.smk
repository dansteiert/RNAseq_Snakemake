rule create_links:
    input:
        metadata = metadata_file,
    output:
        r1 = expand(os.path.join(work_dir, symlink_dir, "{sample}_R1.fastq.gz"), sample=sample),
        r2 = expand(os.path.join(work_dir, symlink_dir, "{sample}_R2.fastq.gz"), sample=sample),
    log:
        os.path.join(work_dir, log_dir, "symlink.log")
    resources:
        mem_mb=1000,
        runtime=20,
        nodes=1,
    threads: 1
    run:
        import pandas as pd
        import os
        df_metadata = pd.read_csv(input.metadata)
        os.makedirs(os.path.join(work_dir, symlink_dir), exist_ok=True)
        for index, row in df_metadata.iterrows():
            # Fastq File
            fastq_file = Path(row["Absolute_Path"])
            output_file_Read = os.path.join(work_dir, symlink_dir, f'{row["Sample_Name"]}_{row["Read"]}.fastq.gz')
            if os.path.exists(output_file_Read):
                os.remove(output_file_Read)
            os.symlink(fastq_file, output_file_Read)
        print(os.listdir(os.path.join(work_dir, symlink_dir)))