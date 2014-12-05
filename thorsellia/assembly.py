import sh
import os
from pandas import DataFrame
from numpy import matrix
from numpy import append
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
import imp
from Bio import SeqIO
from pandas import MultiIndex
from Bio import Entrez
import threading
import time
import sys
from collections import Counter
import json
from thorsellia.clustering import Clustering
from thorsellia.blast import blast
from thorsellia.blast import make_db
from thorsellia.blast import blast_dbs
from thorsellia.parameters import *
from tqdm import tqdm

########################################################################################################3


class InterProThread (threading.Thread):
    def __init__(self, assembly):
        threading.Thread.__init__(self)
        self.assembly = assembly
    def run(self):
        print "Starting thread for", self.assembly, "at", time.ctime(time.time())
        self.assembly.clean_proteins()
        self.assembly.annotate_interpro()
#        self.assembly.annotation_merger()
#        time.sleep(10)
        print "Exiting thread for", self.assembly, "at", time.ctime(time.time())

        
########################################################################################################3
    
class Assembly(object):

    def __init__(self,name,r1_path,r2_path,spade_out,c_min=None,C_max=None,l_min=None, type_scaffold =None, genus = "Thorsellia", species = "anophelis", kingdom = "Bacteria", strain = None, out_path = None):
        self.name = name
        self.genus = genus
        self.species = species
        self.kingdom = kingdom
        self.r1 = r1_path
        self.r2 = r2_path
        self.spades_out = spade_out + self.name +"/"
        self.c_min = c_min
        self.C_max = C_max
        self.l_min = l_min
        if out_path:
            self.annotation_path = out_path
        else : 
            self.annotation_path = annotation_path
        self.annotation_base = self.annotation_path + self.name + "//" + self.name
        if self.r1:
            self.filtered_scafs = self.spades_out + "scaffolds" + ".filtered.fasta"
        else :
            self.filtered_scafs = spade_out
        self.type_scaffold = type_scaffold
        self.proteins =  self.annotation_base + ".faa"
        self.genes =  self.annotation_base + ".ffn"
        self.clean_built =  self.annotation_base + ".fsa"
        self.ip_base =  self.annotation_base + "_interpro"
        self.prokka_gff = self.annotation_base + ".gff"
        self.interpro_gff = self.annotation_base + "_interpro.gff3"
        self.gff = self.annotation_base + "_full.gff"
        self.scc_path = self.annotation_path + self.name + "//" + "SCCs"
        self.sccs_proteins = self.annotation_base + "_sccs.faa"
        self.sccs_genes = self.annotation_base + "_sccs.ffn"
        self.db_file = self.annotation_base + "_full.db"
        if strain:
            self.strain = strain
        else :
            self.strain = self.name
        
    def assemble(self):
        print "assembling", self.name
        spades_script("-o",self.spades_out,  "-1", self.r1, "-2", self.r2,  "-t", threads, "--careful") 

    def filter(self):
        print "filtering", self.name
        filterer = imp.load_source('spades_filter', '/home/moritz/repos/Pyscratches/20140530_spades_filters/spades_fasta_stats.py')
        return filterer.filter_assembly(self.spades_out + "scaffolds.fasta" , self.spades_out + "scaffolds", auto=False, c_min=self.c_min, C_max = self.C_max, l_min = self.l_min)

    def dot_plot(self,path = None):
        if not path: path = spades_out + "dot_plot"
        print "dot-plot", self.name
        sh.nucmer("-D", 4 , "-d", 0.08, "-g", 500, "-l", 8, "-c", 30, "-b", 600,"-p", path,  self.clean_built, self.type_scaffold)
        sh.mummerplot("--color", "--layout", "--postscript", "--large", "-R",  self.clean_built, "-Q", self.type_scaffold, "-p", path,  path + ".delta")

    def annotate_prokka(self):
        print "annotate", self.name
        sh.prokka( "--outdir", self.annotation_path + self.name, "--prefix",  self.name, "--compliant", "--genus", self.genus, "--species", self.species, "--strain", self.strain, "--kingdom", self.kingdom, "--locustag", self.name, "--force", "--cpus", threads,   self.filtered_scafs)       

    def annotate_interpro(self):
        print "Running interpro scan on", self.name
        interpro_script("--output-file-base", self.ip_base, "-dp", "-f", "GFF3", "--goterms", "-i", self.proteins, "-pa")
        print("--output-file-base", self.ip_base, "-dp", "-f", "GFF3", "--goterms", "-i", self.proteins, "-pa")

    def annotation_merger(self):
        print "merging gffs on", self.name
        GFFmerger = imp.load_source('GFFmerger', '/home/moritz/repos/Pyscratches/20140623_gffmerger/gffmerger.py')
        annotation_db = GFFmerger.GFFmerger(self.prokka_gff,self.interpro_gff, self.db_file)
        annotation_db.write_db(self.gff)

    def clean_proteins(self):
        print "Cleaning protein file of", self.name
        sh.sed("-i", 's/*$//', self.proteins)

    def find_SCCs(self):
        scc_dbs = {d : blast_dbs[d] for d in blast_dbs if "cog" in blast_dbs[d]['sets'] and "SCG" in blast_dbs[d]['sets']}
        db_max_size = max([os.stat(scc_dbs[db]['path']).st_size for db in scc_dbs])
        search_space_size = db_max_size + os.stat(self.genes).st_size

        results = {}
        for db in tqdm(scc_dbs) :
            blast_out = blast(self.genes, self.scc_path, scc_dbs[db]['path'], "blastn", clean = True, post_process = False, evalue = 0.01, search_space = search_space_size)
            results[db] = blast_out.iloc[list(blast_out['length'] == blast_out['length'].max())].iloc[0]

        dico = {v.name[0] : k for k,v in results.iteritems()}
        with open(self.proteins) as files:
            proteins = [s for s in SeqIO.parse(files,"fasta")]

        proteins = [p for p in proteins if p.id in dico.keys()]
        for p in proteins:
            p.id = dico[p.id]
            p.description = ""

        with open(self.sccs_genes, "w") as files:
            SeqIO.write(proteins,files,"fasta")

        with open(self.genes) as files:
            genes = [s for s in SeqIO.parse(files,"fasta")]
            
        genes = [p for p in genes if p.id in dico.keys()]
        for p in genes:
            p.id = dico[p.id]
            p.description = ""

        with open(self.sccs_genes, "w") as files:
            SeqIO.write(genes,files,"fasta")
            
        
