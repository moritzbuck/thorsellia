import sh
import os
from Bio import SeqIO
from thorsellia.parameters import *
from Bio import Seq

seed = 23

def mask_seq(seq,mask):
    chars = str(seq.seq)
    chars = [c for i,c in enumerate(chars) if mask[i]]
    out = SeqIO.SeqRecord(Seq.Seq("".join(chars)))
    out.id = seq.id 
    return out

class SixteenS(object):

    # got from arb a file with the 50 most similar 16 to otu_4 the cat of the fasta is the input
    
    def __init__(self, name, path , raw):
        self.name = name
        with open(raw, "r") as ff:
            seqs = [s for s in SeqIO.parse(ff,"fasta")]
        self.seqs = [s for s in seqs if  len(s) > 395 ]
        unique_seqs = list(set([str(s.seq) for s in seqs]))
    
        cleaned_seqs = []
        for s in self.seqs:
                if str(s.seq) in unique_seqs:
                    cleaned_seqs.append(s)
                    unique_seqs.remove(str(s.seq))
        self.seqs = cleaned_seqs
        
        with open(path + name + ".filtered.fasta", "w") as ff:
                SeqIO.write(self.seqs,ff,"fasta")
        self.raw_path = path + name + ".filtered.fasta"
        self.root_path = path
        self.root =  self.root_path + self.name
        self.seed = 23
        
        
    def align(self):
        print "muscle it (e.g. align)"
        
        sh.muscle("-in", self.raw_path,"-out", self.root + ".aligned.fasta","-maxiters", 1 ,"-diags")

        print "change names"
        with open(self.root + ".aligned.fasta","r") as file:
            self.seqs = [s for s in SeqIO.parse(file,"fasta") ]
            for s in self.seqs:
                s.description = ""


        with open(self.root + ".aligned.fasta","w") as outp:
            SeqIO.write(self.seqs,outp,"fasta")

    def block(self):
        with open(self.root + ".aligned.fasta","r") as file:
            seqs = [s for s in SeqIO.parse(file,"fasta")]
        otu_4 = [s for s in seqs if s.id =="OTU_4"][0]
        mask = [c != "-" for c in str(otu_4.seq)]
        seqs_masked = [mask_seq(s,mask) for s in seqs]
        seq_covered = [ s  for s in seqs_masked if float(len(s)-str(s.seq).count("-"))/len(s) > 0.1 ]
        seqs = seq_covered
        for s in seqs:
            s.description = ""

        unique_seqs = list(set([str(s.seq) for s in seqs]))
    
        cleaned_seqs = []
        for s in seqs:
            if str(s.seq) in unique_seqs:
                cleaned_seqs.append(s)
                unique_seqs.remove(str(s.seq))

        seqs = cleaned_seqs
        papers_w = [s for s in seqs if "gi" in s.id][0]
        
        mask = [c != "-" for c in str(papers_w.seq)]
        seqs_masked = [mask_seq(s,mask) for s in seqs]
        seq_covered = [ s  for s in seqs_masked if float(len(s)-str(s.seq).count("-"))/len(s) == 1 ]
        seqs = seq_covered
        for s in seqs:
            s.description = ""

        unique_seqs = list(set([str(s.seq) for s in seqs]))
    
        cleaned_seqs = []
        for s in seqs:
            if str(s.seq) in unique_seqs:
                cleaned_seqs.append(s)
                unique_seqs.remove(str(s.seq))
        
        with open(self.root + ".aligned.cleaned.fasta","w") as outp:
            SeqIO.write(cleaned_seqs,outp,"fasta")

        print "make a block"
        try:
            sh.Gblocks(self.root + ".aligned.cleaned.fasta", "-tD")
        except:
            pass
        
        os.remove(self.root + ".aligned.cleaned.fasta")
        if os.path.exists(self.root + ".aligned.cleaned.fasta" + "-gb.htm"):
            os.remove(self.root + ".aligned.cleaned.fasta" + "-gb.htm")
            os.rename(self.root + ".aligned.cleaned.fasta" + "-gb", self.root + ".aligned.fasta")
            print "remove identical seqs"

            with open(self.root + ".aligned.fasta","r") as file:
                seqs = [s for s in SeqIO.parse(file,"fasta")]

            unique_seqs = list(set([str(s.seq) for s in seqs]))
    
            print len(unique_seqs)
            cleaned_seqs = []
            for s in seqs:
                if str(s.seq) in unique_seqs:
                    cleaned_seqs.append(s)
                    unique_seqs.remove(str(s.seq))
    
            print len(unique_seqs)

            with open(self.root + ".aligned.fasta","w") as outp:
                SeqIO.write(cleaned_seqs,outp,"fasta")

            return 1
        else : 
            return 0

        
    def tree_construction(self):
        print "build a tree"
        if os.path.exists(self.root + "RAxML/" ):
            sh.rm("-r", self.root + "RAxML/")
        os.makedirs(self.root + "RAxML/")

        sh.raxmlHPC_PTHREADS_AVX("-w", self.root + "RAxML/", "-T", threads-2, "-m", "GTRGAMMA", "-p", self.seed, "-#", 20, "-s", self.root + ".aligned.fasta", "-n", "T13") 
        print "boostrap dat tree"
        sh.raxmlHPC_PTHREADS_AVX("-w", self.root + "RAxML/", "-T", threads-2, "-m", "GTRGAMMA", "-p", self.seed, "-b", self.seed, "-#", 100, "-s", self.root + ".aligned.fasta", "-n", "T14")
        print "combine"
        sh.raxmlHPC_AVX("-m", "GTRGAMMA", "-w", self.root + "RAxML/", "-p", self.seed, "-f", "b", "-t", self.root + "RAxML/"+"RAxML_bestTree.T13", "-z",self.root + "RAxML/"+ "RAxML_bootstrap.T14", "-n", "T15")
        print "clean up"
#        os.remove(self.root + "_branches_labeled.tree")
#        os.remove(self.root + "_nodes_labeled.tree")
#        sh.ln("-s",  self.root + "RAxML/RAxML_bipartitionsBranchLabels.T15", self.root +"_branches_labeled.tree")
#        sh.ln("-s",  self.root + "RAxML/RAxML_bipartitions.T15", self.root +"_nodes_labeled.tree")

