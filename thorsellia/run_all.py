#!/usr/bin/env python
#SBATCH -A b2013086
#SBATCH  -p node
#SBATCH  -n 16
#SBATCH  --mail-user murumbii@gmail.com
#SBATCH  -t 1-0:00:00 

import thorsellia
from thorsellia import Assembly



bowtie2_out_path = type_correction_path + "bowtie/"
bowtie2_idx_name = "thorsellia_type"
bowtie2_out_sam = bowtie2_out_path + bowtie2_idx_name + "_mapped.sam"
bowtie2_out_bam = bowtie2_out_path + bowtie2_idx_name + "_mapped.bam"

os.chdir(bowtie2_out_path)
sh.bowtie2_build(type_scaffold, bowtie2_idx_name)
sh.bowtie2("-x", bowtie2_out_path + bowtie2_idx_name, "-1", type_reads_r1, "-2", type_reads_r2, "-S", bowtie2_out_sam, "-N", 1, "--local", "-p", threads )
sh.samtools("view","-bS", bowtie2_out_sam, "-o" , bowtie2_out_bam, "-@", threads)
sh.samtools("sort", bowtie2_out_bam , bowtie2_out_bam.split(".")[0])
sh.samtools("index", bowtie2_out_bam)

depths = array([int(line.split("\t")[2]) for line in sh.samtools("depth", bowtie2_out_bam, _iter = True )])
depths_hist = histogram(depths, bins=100)

pdf = PdfPages('coverage_distrib.pdf')


"""
Using SEQuel (http://bix.ucsd.edu/SEQuel) to correct the type scaffold with the new reads
Needs BWA in path
Versions used:
SEQuel v1.0.2
BWA 0.7.8
"""

SEQuel_output = work_root + "type_strain_correction"
SEQuel_path = "/home/moritz/src/SEQuel/"
SEQuel_prep_path = SEQuel_path + "prep.pl"
SEQuel_jar_path = SEQuel_path + "SEQuel.jar"




os.chdir(SEQuel_output)
if os.path.exists("prep"): sh.rm("-r","prep")
SEQuel_prep = sh.Command(SEQuel_prepgen_path)
SEQuel_prep("-r1", type_reads_r1,  "-r2", type_reads_r2, "-t", threads, "-c", type_scaffold)
with open("prep/prep.log") as log:
    log = log.readlines()
    
log = [l for l in log if "insert" in l][0]
insert_length = int(log.split(" ")[-1].split("\n")[0])
sh.java("-Xmx64g", "-jar", SEQuel_jar_path, "-A", "prep", "-i", insert_length, "-p", threads)





#################################################################################################

import illumitag
from illumitag.clustering import Cluster

""" running amplicon pipeline """

#illumitag.projects['louise'].run_pools_slurm()


louise_samples = []

for p in illumitag.projects['louise']: louise_samples += p.samples

louise_samples = [l for l in louise_samples if l.used]
louise_clusters = Cluster(louise_samples,"louise")

louise_clusters.process_samples()
louise_clusters.combine_reads()
louise_clusters.otu_uparse.run()
louise_clusters.otu_uparse.taxonomy_silva.assign()
louise_clusters.otu_uparse.taxonomy_silva.make_otu_table()
louise_clusters.otu_uparse.taxonomy_silva.make_otu_table_norm()
louise_clusters.otu_uparse.taxonomy_silva.make_plots()
louise_clusters.otu_uparse.taxonomy_silva.stats.nmds.run()
louise_clusters.otu_uparse.taxonomy_silva.comp_phyla.make_taxa_table()
louise_clusters.otu_uparse.taxonomy_silva.comp_phyla.make_plots()
louise_clusters.otu_uparse.taxonomy_silva.comp_phyla.stats.nmds.run()
louise_clusters.otu_uparse.taxonomy_silva.comp_tips.make_taxa_table()
louise_clusters.otu_uparse.taxonomy_silva.comp_tips.make_plots()
louise_clusters.otu_uparse.taxonomy_silva.comp_tips.stats.nmds.run()
louise_clusters.otu_uparse.taxonomy_silva.comp_class.make_taxa_table()
louise_clusters.otu_uparse.taxonomy_silva.comp_class.make_plots()
louise_clusters.otu_uparse.taxonomy_silva.comp_class.stats.nmds.run()
louise_clusters.otu_uparse.taxonomy_silva.comp_order.make_taxa_table()
louise_clusters.otu_uparse.taxonomy_silva.comp_order.make_plots()
louise_clusters.otu_uparse.taxonomy_silva.comp_order.stats.nmds.run()

louise_clusters_99 = Cluster(louise_samples,"louise_99")

louise_clusters_99.process_samples()
louise_clusters_99.combine_reads()
louise_clusters_99.otu_uparse.run(threshold=1)
louise_clusters_99.otu_uparse.taxonomy_silva.assign()
louise_clusters_99.otu_uparse.taxonomy_silva.make_otu_table()
louise_clusters_99.otu_uparse.taxonomy_silva.make_otu_table_norm()
louise_clusters_99.otu_uparse.taxonomy_silva.make_plots()
louise_clusters_99.otu_uparse.taxonomy_silva.stats.nmds.run()
louise_clusters_99.otu_uparse.taxonomy_silva.comp_phyla.make_taxa_table()
louise_clusters_99.otu_uparse.taxonomy_silva.comp_phyla.make_plots()
louise_clusters_99.otu_uparse.taxonomy_silva.comp_phyla.stats.nmds.run()
louise_clusters_99.otu_uparse.taxonomy_silva.comp_tips.make_taxa_table()
louise_clusters_99.otu_uparse.taxonomy_silva.comp_tips.make_plots()
louise_clusters_99.otu_uparse.taxonomy_silva.comp_tips.stats.nmds.run()
louise_clusters_99.otu_uparse.taxonomy_silva.comp_class.make_taxa_table()
louise_clusters_99.otu_uparse.taxonomy_silva.comp_class.make_plots()
louise_clusters_99.otu_uparse.taxonomy_silva.comp_class.stats.nmds.run()
louise_clusters_99.otu_uparse.taxonomy_silva.comp_order.make_taxa_table()
louise_clusters_99.otu_uparse.taxonomy_silva.comp_order.make_plots()
louise_clusters_99.otu_uparse.taxonomy_silva.comp_order.stats.nmds.run()
