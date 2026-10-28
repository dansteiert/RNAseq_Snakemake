rule generate_count_matrix:
    input:
        alignment = expand(os.path.join(work_dir, alignment_dir, "{sample}_ReadsPerGene.out.tab"), sample=sample)
    output:
        count_Matrix = os.path.join(work_dir, count_matrix_dir, "CountMatrix.csv"),
        count_Matrix_non_zero = os.path.join(work_dir, count_matrix_dir, "CountMatrix_non_zero_count.csv"),
        count_Matrix_protein_coding_non_zero = os.path.join(work_dir, count_matrix_dir, "CountMatrix_non_zero_count_protein_coding_only.csv")
    params:
        work_dir = work_dir,
        count_matrix_dir = count_matrix_dir,
        sequencing_strand = sequencing_strand,
        samples = sample  # Add samples to params
    log:
        os.path.join(work_dir, log_dir, "gen_count_matrix.log")
    resources:
        mem_mb=4000,
        runtime=60,
        nodes=1
    threads: 1
    run:
        import pandas as pd
        import os

        # Determine Strandedness from Input
        use_col_translation_dict = {"unstranded": 1, "forward_stranded": 2, "reverse_stranded": 3}
        quant_col = use_col_translation_dict.get(params.sequencing_strand, "FAILED")
        if quant_col == "FAILED":
            print(params.sequencing_strand, "is not in", use_col_translation_dict.keys())
            return

        # Create list of dataframes with correct sample IDs
        df_collection = []
        for file, sample_id in zip(input.alignment, params.samples):
            df = pd.read_csv(file,
                           sep="\t",
                           skiprows=4,
                           index_col=[0],
                           header=None,
                           usecols=[0, quant_col],
                           names=["index", sample_id])
            df_collection.append(df)

        # Add Gene Infos to the Table
        first_file = input.alignment[0]
        geneinfo_path = first_file.replace("_ReadsPerGene.out.tab", "__STARgenome/geneInfo.tab")
        df_gene_info_cols = pd.read_csv(geneinfo_path,
                                      sep="\t",
                                      skiprows=1,
                                      index_col=[0],
                                      names=["Ensembl_ID", "Symbol", "Type"])
        
        # Append gene info and concatenate all dataframes
        df_collection.append(df_gene_info_cols)
        df_out = pd.concat([df_collection[-1], *df_collection[:-1]],  # set the last as the first df such that the identifiers are in the front
                          axis=1,
                          join="outer")

        # Write out different Configuration Files
        # Full Table
        df_out.to_csv(output.count_Matrix)

        # Non Zero Counts
        count_cols = [col for col in df_out.columns if col not in {"Symbol", "Type"}]
        df_out = df_out[df_out[count_cols].sum(axis=1) > 0]
        df_out.to_csv(output.count_Matrix_non_zero)

        # Protein Coding Only
        df_out = df_out[df_out["Type"] == "protein_coding"]
        df_out.to_csv(output.count_Matrix_protein_coding_non_zero) 