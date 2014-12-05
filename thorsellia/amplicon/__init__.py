import illumitag
from illumitag.clustering import Cluster

class Amplicon(object):

    def __init__(self,name = 'louise' )
        self.name = name
        self.pools = illumitag.projects[self.name]
        
    def run_pools(self):
        self.pools.run_pools_slurm()

    def clustering(self, percent=1):
        c_id = self.name + "_" + str(percent)
        if(percent == 1) c_id = self.name

        samples = []
    
        for p in self.pools: samples += p.samples

        samples = [l for l in samples if l.used]
        clusters = Cluster(samples,c_id)
    
        clusters.process_samples()
        clusters.combine_reads()
        clusters.otu_uparse.run()
        clusters.otu_uparse.taxonomy_silva.assign()
        clusters.otu_uparse.taxonomy_silva.make_otu_table()
        clusters.otu_uparse.taxonomy_silva.make_otu_table_norm()
        clusters.otu_uparse.taxonomy_silva.make_plots()
        clusters.otu_uparse.taxonomy_silva.stats.nmds.run()
        clusters.otu_uparse.taxonomy_silva.comp_phyla.make_taxa_table()
        clusters.otu_uparse.taxonomy_silva.comp_phyla.make_plots()
        clusters.otu_uparse.taxonomy_silva.comp_phyla.stats.nmds.run()
        clusters.otu_uparse.taxonomy_silva.comp_tips.make_taxa_table()
        clusters.otu_uparse.taxonomy_silva.comp_tips.make_plots()
        clusters.otu_uparse.taxonomy_silva.comp_tips.stats.nmds.run()
        clusters.otu_uparse.taxonomy_silva.comp_class.make_taxa_table()
        clusters.otu_uparse.taxonomy_silva.comp_class.make_plots()
        clusters.otu_uparse.taxonomy_silva.comp_class.stats.nmds.run()
        clusters.otu_uparse.taxonomy_silva.comp_order.make_taxa_table()
        clusters.otu_uparse.taxonomy_silva.comp_order.make_plots()
        clusters.otu_uparse.taxonomy_silva.comp_order.stats.nmds.run()

