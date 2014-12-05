
from thorsellia.blast import blast
from thorsellia.blast import make_db
from thorsellia.blast import blast_head

from thorsellia.parameters import *
#
import os
import sh
import json
from Bio import SeqIO
from collections import Counter
from tqdm import tqdm
import sys
import shutil
from pandas import DataFrame

class GeneCluster(object):

    def __repr__(self): return '<%s object %s with %i genes from %i genomes>' % (self.__class__.__name__, self.annotation, len(self.genes), len(self.genomes))

    def __init__(self, clustering, genes ):

        self.clustering = clustering

        if type(genes) == dict :
            self.from_dict(genes)
        else:
            self.from_list(genes)
           
    def to_dict(self):
        return {u'annot_fraction': self.annot_fraction, u'annotation': self.annotation,  u'genes': self.genome_2_gene_map,  u'mapping': self.mapping}

    def from_dict(self,imp):
        self.genes = imp['mapping'].keys()
        self.genomes = imp['genes'].keys()
        self.genome_2_gene_map =  imp['genes']
        
        self.annotation = imp['annotation']
        self.annot_fraction = imp['annot_fraction']
        self.mapping = imp['mapping']
        
    def from_list(self, genes):
        self.genes = genes
        self.genomes = list(set(["_".join(g.split('_')[:-1]) for g in genes]))
        self.genome_2_gene_map =  {go : [ge for ge in genes if go in ge] for go in self.genomes}
        
        sub_dict = {g : self.clustering.id2name_map[g] for g in self.genes}
        name_counts = Counter(sub_dict.values())
        total = sum([name_counts[z] for z in name_counts])
        annot_frac = float(name_counts.most_common()[0][1])/float(total)

        self.annotation = name_counts.most_common(1)[0][0]
        self.annot_fraction = annot_frac
        self.mapping = sub_dict

    def to_sequences(self):
        with open(self.clustering.db, "r") as db: seqs = [s for s in SeqIO.parse(db, "fasta") if s.id in self.genes]
            
        return seqs

    def calc_checksum(self, s):
        return str(sum(ord(c)*i for i,c in enumerate(s)))

    def align(self,output_file):
        with open("temp.faa","w") as unalign:
            temp_seqs = self.to_sequences()
            for s in temp_seqs:
                s.id = self.clustering.assembly_id("_".join(s.id.split("_")[:-1]) + "_".join(s.id.split("_")[:-1]))+"_"+s.id.split("_")[-1]
                s.description = ""
            SeqIO.write(temp_seqs, unalign, "fasta")
        sh.muscle("-in", "temp.faa","-out", "temp_aligned.faa")
        os.remove("temp.faa")
        try:
            if self.clustering.seq_type == "genes":
                sh.Gblocks("temp_aligned.faa", "-tD")
            else:
                sh.Gblocks("temp_aligned.faa", "-tP")
        except:
            pass
        
        os.remove("temp_aligned.faa")
        if os.path.exists("temp_aligned.faa-gb.htm"):
            os.remove("temp_aligned.faa-gb.htm")
            shutil.move("temp_aligned.faa-gb", output_file)
            return 1
        else : 
            return 0

   


class Clustering(object):

    def __repr__(self): return '<%s object "%s with %i clusters">' % (self.__class__.__name__, self.name, len(self))
    def __iter__(self): return iter(self.clusters)
    def __len__(self): return len(self.clusters)
    def __getitem__(self, key): return self.clusters[key]

        
    def __init__(self, assemblies, path, name, seq_type="proteins",folder = None):

        if folder:
            from thorsellia.assembly import Assembly
            genomes = os.listdir(folder)
            self.assemblies = [Assembly(g,None, None, folder + g + "/" + os.listdir(folder + g)[0], genus = g.split("_")[0], species =  g.split("_")[1], out_path = path, kingdom = "Archea")  for g in genomes]
            for g in tqdm(self.assemblies):
                if not os.path.exists(g.prokka_gff):
                    print "annotating" , g.name
                    g.annotate_prokka()
        else:
            self.assemblies = assemblies
        self.seq_type = seq_type
        self.path = path

        if not os.path.exists(path):
            os.makedirs(path)

        self.seed = 23
        self.name = name
        self.base = self.path + self.name


        self.db = self.base +".fasta"
        self.raw_blast = self.base +".raw.blast"
        self.blast = self.base +".blast"
        self.abc = self.base +".abc"
        self.mci = self.base +".mci"
        self.tab = self.base +".tab"
        self.raw_clusters = self.base + ".cluster"
        self.processed_clusters = self.base + ".json"
        self.align_path = self.base + "_align/"
        self.scc_align_path = self.base + "_scc_align/"
        
        if os.path.exists(self.processed_clusters):
            with open(self.processed_clusters) as file:
                self.clusters= [GeneCluster(self,c) for c in json.load(file)]
                
        if os.path.exists(self.db):
            with open(self.db,"r") as fas:
                self.id2name_map = {s.description.split(" ")[0] : " ".join(s.description.split(" ")[1:]) for s in SeqIO.parse(fas, "fasta")}

    def assembly_id(self, a):
        return str(sum(ord(c)*i for i,c in enumerate(a)))

                
    def do_self_blasting(self):
        make_db(self.assemblies, self.db, seq_type = self.seq_type)
        blast(self.db, self.base, db = self.db , alg = "blastp" if self.seq_type == "proteins" else "blastn", clean=False, post_process=False, outfmt=6, evalue = 10, word_size = 2 if self.seq_type == "proteins" else 9)

    def do_blast_filter(self):
        with open(self.raw_blast) as inpt:
            if(self.seq_type == "proteins"):
                sh.filtersearchio5("-qcoverage", 70, "-identity", 40, "-format", 8, _in = inpt, _out = self.blast)
            else:
                sh.filtersearchio5("-qcoverage", 50, "-identity", 30, "-format", 8, _in = inpt, _out = self.blast)

    def do_blast_parsing(self, qcov = 0.50, identity = 30):
        raw_blast = DataFrame.from_csv(self.raw_blast, header=-1,index_col=[0,1],sep="\t")
        raw_blast.columns = blast_head[2:12]
        raw_blast.index.names = blast_head[0:2]

        with open(self.db) as file:
            dicti = {s.id : len(s) for s in SeqIO.parse(file,"fasta")}
        covs = (raw_blast['subjectend']-raw_blast['subjectstart']+1)/ [float(dicti[n]) for n in raw_blast.index.get_level_values(0)]
        filtered_blast = raw_blast[covs > qcov]
        # -raw_blast['mismatches']- raw_blast['gapopenings']
        #        filtered_blast = raw_blast[[float(r[1][5]-r[1][4]+1-r[1][2]-r[1][3]) / float(dicti[r[1].name[0]]) > qcov  for r in raw_blast.iterrows()]]
        filtered_blast = filtered_blast[filtered_blast['identity'] > 30]

 #       1/0
        filtered_blast.to_csv(self.blast,sep="\t", header=False)

    def do_mcl(self, inflation = 1.5, resource = 4):
        sh.cut("-f", "1,2,11", self.blast, _out = self.abc)
        sh.mcxload("-abc", self.abc, "--stream-mirror", "--stream-neg-log10", "-stream-tf", "ceil(200)", "-o",  self.mci, "-write-tab"  , self.tab)
        sh.mcl(self.mci, "-I", inflation , "-resource", resource, "-use-tab",  self.tab, "-o",  self.raw_clusters)
    
    def single_copy_clusters(self):
        return [c for c in self.clusters if len(c.genes) == len(self.assemblies) and len(c.genomes)==len(self.assemblies)]

    def almost_single_copy_clusters(self):
        return [c for c in self.clusters if len(c.genes) == (len(self.assemblies)-1) and (len(c.genomes)==len(self.assemblies)-1)]
    

    def post_process(self):
        with open(self.raw_clusters) as c_file:
            clusters=[l[:-1].split("\t") for l in c_file.readlines()]
            clusters = [[g.split("(")[0] for g in v ] for v in clusters]
        self.clusters = []

        with open(self.db,"r") as fas:
            self.id2name_map = {s.description.split(" ")[0] : " ".join(s.description.split(" ")[1:]) for s in SeqIO.parse(fas, "fasta")}
                    
        for i,c in enumerate(clusters):
            self.clusters.append(GeneCluster(self,c))
            
        with open(self.processed_clusters, 'w') as outfile:
            json.dump([c.to_dict() for c in self.clusters], outfile,  indent=4, sort_keys=True)

    def run(self):
        print "Doing MCL clustering"

        print "blasting"
        self.do_self_blasting()
        
        print "filter blast results"
#        self.do_blast_filter()
        if(self.seq_type == "proteins"):
            self.do_blast_parsing( qcov = 0.5, identity = 30)
        else :
            self.do_blast_parsing()
                              
        print "do the MCL-pipe"
        self.do_mcl()

        print "post-process"
        self.post_process()
        
    def align_all(self):
        print "aligning the hell out of it"
        if os.path.exists(self.align_path):
             shutil.rmtree(self.align_path)
        os.makedirs(self.align_path)

        for i,clust in tqdm(enumerate(self.single_copy_clusters())):
            clust.align(self.align_path + str(i) + ".fasta")

        
    def cat_align(self):
        print "CONCATENATE!!!11"
        cat_seqs = None
        for f in os.listdir(self.align_path):
            if "fasta" in f:
                with open(self.align_path + f,"r") as file:
                    seqs = [s for s in SeqIO.parse(file,"fasta")]
                    if not cat_seqs:
                        cat_seqs = seqs
                        order = ["_".join(s.id.split("_")[:-1]) for s in cat_seqs]
                        for i,s in enumerate(cat_seqs):
                            s.id = order[i]
                    else :
                        seqs_sorted = [[z for z in seqs if o in z.id][0] for o in order]
                        cat_seqs = [ s + seqs_sorted[i] for i,s in enumerate(cat_seqs) ]

        a2id = {a.name : self.assembly_id(a.name + a.name)  for a in self.assemblies}
        id2a = {a2id[a] : a for a in a2id}
        for i,s in enumerate(cat_seqs):
            s.id = id2a[order[i]]
            s.description = "composite of " + str(len([f for f in os.listdir(self.align_path) if "fasta" in f])) + " single copy genes-clusters"
        with open(self.base + "_cat_align.fasta","w") as outp:
            SeqIO.write(cat_seqs,outp,"fasta")
        return cat_seqs
        
    def tree_construction(self,root = None, sccs = False):
        print "build a tree"
        if os.path.exists(self.base + "RAxML/" ):
            sh.rm("-r", self.base + "RAxML/")
        os.makedirs(self.base + "RAxML/")

        if self.seq_type == "proteins" :
            model = "PROTGAMMAAUTO"
        else:
            model = "GTRGAMMA"

        alignment = self.base + "_scc_cat_align.fasta" if sccs else self.base + "_cat_align.fasta"
        
        sh.raxmlHPC_PTHREADS_AVX("-w", self.base + "RAxML/", "-T", threads-2, "-m", model, "-p", self.seed, "-#", 20, "-s", alignment, "-n", "T13")#, "-o", root) 
        print "boostrap dat tree"
        sh.raxmlHPC_PTHREADS_AVX("-w", self.base + "RAxML/", "-T", threads-2, "-m", model, "-p", self.seed, "-b", self.seed, "-#", 100, "-s", alignment, "-n", "T14")#, "-o", root)
        print "combine"
        sh.raxmlHPC_AVX("-m", "GTRCAT", "-w", self.base + "RAxML/", "-p", self.seed, "-f", "b", "-t", self.base + "RAxML/"+"RAxML_bestTree.T13", "-z",self.base + "RAxML/"+ "RAxML_bootstrap.T14", "-n", "T15")#, "-o", root)
        print "clean up"
        if os.path.exists(self.base + "_branches_labeled.tree"):
            os.remove(self.base + "_branches_labeled.tree")
            os.remove(self.base + "_nodes_labeled.tree")
        sh.ln("-s",  self.base + "RAxML/RAxML_bipartitionsBranchLabels.T15", self.base +"_branches_labeled.tree")
        sh.ln("-s",  self.base + "RAxML/RAxML_bipartitions.T15", self.base +"_nodes_labeled.tree")

        
    def rm_genome(self, name):
        self.assemblies = [a for a in self.assemblies if name not in a.name]
        for c in self:                                                             
            c.genomes = [g for g in  c.genomes if name not in g]
        for c in self:                                                             
            c.genes = [g for g in  c.genes if name not in g]


    def keep_genomes(self, genomes):
        self.assemblies = [a for a in self.assemblies if a.name in genomes]
        for c in self:                                                             
            c.genomes = [g for g in  c.genomes if g in genomes ]
        for c in self:                                                             
            c.genes = [g for g in  c.genes if len([o for o in genomes if g.count(o) ==1]) > 0 ]


    def cooccurence_matrix(self):
        matrix = DataFrame(data=0, index = [a.name for a in self.assemblies], columns=[a.name for a in self.assemblies])
        for a in self.assemblies:
            for b in self.assemblies:
                count = 0
                for c in self:
                    if a.name in c.genomes and b.name in c.genomes:
                        count += 1
                matrix[a.name][b.name] = count
                
        return matrix

    def sccs_cat_align(self):
        all_sccs = {}
        for a in self.assemblies:
            with open(a.sccs_genes) as fasta:
                all_sccs[a.name] = [s for s in SeqIO.parse(fasta,"fasta")]
        cogs = list(set(sum([[s.id for s in seqs] for seqs in all_sccs.values()],[])))
        cogs.sort()

        a2id = {a.name : self.assembly_id(a.name + a.name)  for a in self.assemblies}
        id2a = {a2id[a] : a for a in a2id}
        
        for c in cogs:
            c_seqs = []
            for a,seqs in all_sccs.iteritems():
                seq = [s for s in seqs if s.id == c]
                if len(seq) == 1:
                    seq = seq[0]
                else:
                    continue
                seq.id = a2id[a]
                seq.description = ""
                c_seqs.append(seq)
            if len(c_seqs) == len(self.assemblies):
                with open("temp.faa","w") as fasta:
                    SeqIO.write(c_seqs,fasta,"fasta")
                print "aligning genes for",c
                sh.muscle("-in", "temp.faa","-out", "temp_aligned.faa")
                os.remove("temp.faa")
                try:
                    sh.Gblocks("temp_aligned.faa", "-tD")
                except:
                    pass
        
                os.remove("temp_aligned.faa")
                if os.path.exists("temp_aligned.faa-gb.htm"):
                    os.remove("temp_aligned.faa-gb.htm")
                    if not os.path.exists(self.scc_align_path):
                        os.makedirs(self.scc_align_path)
                    shutil.move("temp_aligned.faa-gb", self.scc_align_path + c + ".ffn")


        print "CONCATENATE!!!11"
        cat_seqs = None
        for f in os.listdir(self.scc_align_path):
            if "ffn" in f:
                with open(self.scc_align_path + f,"r") as file:
                    seqs = [s for s in SeqIO.parse(file,"fasta")]
                    if not cat_seqs:
                        cat_seqs = seqs
                        order =  [s.id for s in cat_seqs]
                        for i,s in enumerate(cat_seqs):
                            s.id = order[i]
                    else :
                        seqs_sorted = [[z for z in seqs if o in z.id][0] for o in order]
                        cat_seqs = [ s + seqs_sorted[i] for i,s in enumerate(cat_seqs) ]
            for i,s in enumerate(cat_seqs):
                s.id = id2a[order[i]]
                s.description = "composite of " + str(len([f for f in os.listdir(self.scc_align_path) if "ffn" in f])) + " single copy cogs"
            with open(self.base + "_scc_cat_align.fasta","w") as outp:
                SeqIO.write(cat_seqs,outp,"fasta")
                
        return all_sccs
