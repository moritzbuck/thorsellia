from sh import Command

threads = 16

spades_script = Command("spades.py")
#interpro_script = Command("/home/moritz/src/interproscan-5.4-47.0/interproscan.sh")

raw_path = "/proj/b2013086/"
work_root = raw_path + "thorsellia/"
type_correction_path = work_root + "type_strain_correction/"
spades_out_path = work_root + "assemblies/"
out_group_path = work_root + "outgroups/"
db_path = work_root + "databases/"
type_scaffold = db_path + "thorsellia_type.fasta"
genes_db = db_path + "all_assemblies_genes.fasta"
genome_db = db_path + "all_assemblies_genomes.fasta"
proteins_db = db_path + "all_assemblies_proteins.fasta"
proteins_of_interest_db = db_path + "proteins_of_interest.fasta"
genes_of_interest_db = db_path + "genomics_of_interest.fasta"
annotation_path = work_root + "annotation/"
gvp_path = annotation_path + "GvPs/"
tmp_folder =  work_root + "temp/"


