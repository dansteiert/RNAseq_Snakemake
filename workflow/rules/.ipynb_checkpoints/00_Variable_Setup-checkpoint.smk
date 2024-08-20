gtf= config["gtf"]
genome= config["genome"]
genome_version = config["genome_version"]
metadata_file = config["metadata_file"] # Sample_Name, Read (R1/R2), Absolute_Path
read_length=config["read_length"]
sample = config["sample"]

work_dir = config["work_dir"]
genome_index_dir = config["genome_index_dir"]

two_pass = config["two_pass"]
sequencing_strand = config["sequencing_strand"] # unstranded, forward, reverse


log_dir = "Log"
symlink_dir = "Symlinks"
star_version = "2.7.11b"
alignment_dir = "Alignment"
count_matrix_dir = "Count_Matrix"