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
        sequencing_strand = sequencing_strand
    log:
        os.path.join(work_dir, log_dir, "gen_count_matrix.log")
    resources:
        mem_mb=4000,
        runtime=60 ,
        nodes=1,
    threads: 1
    run:
        import pandas as pd
        import os
        # Determine Strandedness from Input
        use_col_translation_dict = {"unstranded": 1 , "forward_stranded": 2 , "reverse_stranded": 3}
        quant_col = use_col_translation_dict.get(params.sequencing_strand, "FAILED")
        if quant_col == "FAILED":
            print(use_col, "is not in", use_col_translation_dict.keys())
            return
        #
        df_collection = [pd.read_csv(file,
                                     sep="\t", skiprows=4, index_col=[0], header=None, usecols=[0, quant_col],
                                     names=["index", os.path.split(file)[0].replace("_ReadsPerGene.out.tab", "")])
                         for file in input.alignment]
        # Add Gene Infos to the Table
        df_gene_info_cols = pd.read_csv(os.path.join(input.alignment[0].replace("_ReadsPerGene.out.tab", "__STARgenome"), "geneInfo.tab"),
                                        sep="\t", skiprows=1, index_col=[0],
                                        names=["Ensembl_ID", "Symbol", "Type"])
        df_collection.append(df_gene_info_cols)
        df_out = pd.concat([df_collection[-1], *df_collection[:-1]], # set the last as the first df such that the identifiers are in the front
                           axis=1, join="outer")
        # Write out different Configuration Files
        # Full Table
        df_out.to_csv(os.path.join(params.work_dir, params.count_matrix_dir, "CountMatrix.csv")) 
        # Non Zero Counts
        df_out = df_out[df_out[[i for i in df_out.columns if i not in {"Symbol", "Type"}]].sum(axis=1) > 0]
        df_out.to_csv(os.path.join(params.work_dir, params.count_matrix_dir,
                                   "CountMatrix_non_zero_count.csv"))
        # Protein Coding Only
        df_out = df_out[df_out["Type"]== "protein_coding"]
        df_out.to_csv(os.path.join(params.work_dir, params.count_matrix_dir, "CountMatrix_non_zero_count_protein_coding_only.csv"))    