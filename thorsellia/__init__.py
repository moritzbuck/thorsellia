#!/home/moritz/.pyenv/shims/python

import threading
import time
import sys
import sh
from thorsellia.clustering import Clustering
from thorsellia.blast import blast
from thorsellia.blast import make_db
from thorsellia.parameters import *
from thorsellia.assembly import Assembly
from thorsellia.assembly import InterProThread
from thorsellia.sixteenS import SixteenS

def full_pipe():
    for a in all_assemblies : a.assemble()
    for a in all_assemblies : a.filter()  
    for a in all_assemblies : a.dot_plot()
    make_db(all_assemblies,path = genes_db)
    make_db(all_assemblies,path = genome_db,seq_type = "genomes")
    make_db(all_assemblies,path = proteins_db,seq_type = "proteins")
    blast(proteins_of_interest_db,output_path = annotation_path,  db = proteins_db).to_csv(annotation_path + "Proteins_of_Interest_single_best_hits.csv",clean=False)
    blast(genes_of_interest_db, output_path = annotation_path , db = genome_db, alg = "blastn").to_csv(annotation_path + "Genomics_of_Interest_single_best_hits.csv")
    
    for a in all_assemblies :
        InterProThread(a).start()
    nb_ip_threads = lambda : len([t.__class__.__name__ for t in threading.enumerate() if t.__class__.__name__ == "InterProThread"])
    while nb_ip_threads() > 0:
        print nb_ip_threads(), "InterProScans still running"
        time.sleep(1200)        
    for a in all_assemblies : a.annotation_merger()
        
    for a in all_outgroups : a.annotate_prokka()
        
    clustered_all_nucl.run()
    clustered_all_nucl.align_all()
    clustered_all_nucl.cat_align()
    clustered_all_nucl.tree_construction()    
    return

def run():
#    G_apicola.annotate_prokka()
    clustered_all_nucl.run()
    clustered_all_nucl.align_all()
    clustered_all_nucl.cat_align()
    clustered_all_nucl.tree_construction(root="V_cholerae,P_profundum")    
    return
    

if sh.which("blastp") == None :
    print "no blast avail"
    sys.exit()


Ta   = Assembly("Ta", raw_path   + "INBOX/130411_M00485_0043_000000000-A3EV0/Sample_Ta/Ta_ACTTGA_L001_R1_001.fastq.gz"    , raw_path + "INBOX/130411_M00485_0043_000000000-A3EV0/Sample_Ta/Ta_ACTTGA_L001_R2_001.fastq.gz"    , spades_out_path, c_min=100, C_max=180, l_min=9000, type_scaffold=type_scaffold)
T2_1 = Assembly("T2_1", raw_path + "INBOX/130411_M00485_0043_000000000-A3EV0/Sample_T2_1/T2_1_CAGATC_L001_R1_001.fastq.gz", raw_path + "INBOX/130411_M00485_0043_000000000-A3EV0/Sample_T2_1/T2_1_CAGATC_L001_R2_001.fastq.gz", spades_out_path, c_min=100, C_max=180, l_min=8000, type_scaffold=type_scaffold)
B8   = Assembly("B8", raw_path   + "INBOX/130411_M00485_0043_000000000-A3EV0/Sample_B8/B8_GCCAAT_L001_R1_001.fastq.gz"    , raw_path + "INBOX/130411_M00485_0043_000000000-A3EV0/Sample_B8/B8_GCCAAT_L001_R2_001.fastq.gz"    , spades_out_path, c_min=100, C_max=300, l_min=10000, type_scaffold=type_scaffold)
B3   = Assembly("B3", raw_path   + "INBOX/130411_M00485_0043_000000000-A3EV0/Sample_B3/B3_ACAGTG_L001_R1_001.fastq.gz"    , raw_path + "INBOX/130411_M00485_0043_000000000-A3EV0/Sample_B3/B3_ACAGTG_L001_R2_001.fastq.gz"    , spades_out_path, c_min=90, C_max=160, l_min=6000, type_scaffold=type_scaffold)
s4_1 = Assembly("4_1",  raw_path + "INBOX/130411_M00485_0043_000000000-A3EV0/Sample_4_1/4_1_TGACCA_L001_R1_001.fastq.gz"  , raw_path + "INBOX/130411_M00485_0043_000000000-A3EV0/Sample_4_1/4_1_TGACCA_L001_R2_001.fastq.gz"  , spades_out_path, c_min=10, C_max=100, l_min=6500, type_scaffold=type_scaffold)
s1_1 = Assembly("1_1", raw_path  + "INBOX/130411_M00485_0043_000000000-A3EV0/Sample_1_1/1_1_TTAGGC_L001_R1_001.fastq.gz"  , raw_path + "INBOX/130411_M00485_0043_000000000-A3EV0/Sample_1_1/1_1_TTAGGC_L001_R2_001.fastq.gz"  , spades_out_path, c_min=80, C_max=180, l_min=6000, type_scaffold=type_scaffold)

all_assemblies = [Ta, T2_1, B8, B3, s4_1, s1_1]

E_coli = Assembly("E_coli",  None  , None  , out_group_path + "E_coli.fna" , genus = "Escherichia" , species = "coli" , strain = "ATCC" )
Serratia = Assembly("Serratia",  None  , None  , out_group_path + "Serratia.fna" , genus = "Serratia" , species = "" , strain = "ATCC" )
S_enterica = Assembly("S_enterica",  None  , None  , out_group_path +"S_enterica.fna", genus = "Salmonella" , species = "enterica" , strain = "Typhi Ty2" )
P_multocida = Assembly("P_multocida",  None  , None  , out_group_path + "P_multocida.fna" , genus = "Pasteurella" , species = "multocida" , strain = "3480" )
W_glossinidia = Assembly("W_glossinidia",  None  , None  , out_group_path + "W_glossinidia.fna" , genus = "Wigglesworthia" , species = "glossinidia" , strain = "yale" )
B_aphidicola = Assembly("B_aphidicola",  None  , None  , out_group_path + "B_aphidicola.fna" , genus = "Buchnera" , species = "aphidicola" , strain = "5A" )
A_ferroxidans = Assembly("A_ferroxidans",  None  , None  , out_group_path + "A_ferroxidans.fna" , genus = "Acidithiobacillus" , species = "feroxidans" , strain = "ATCC 23270" )
P_profundum = Assembly("P_profundum",  None  , None  , out_group_path + "P_profundum.fna" , genus = "Photobacterium" , species = "profundum" , strain = "SS9" )
V_cholerae = Assembly("V_cholerae",  None  , None  , out_group_path + "V_cholerae.fna" , genus = "Vibrio" , species = "cholerae" , strain = "MJ" )
G_apicola = Assembly("G_apicola",  None  , None  , out_group_path + "G_apicola.fna" , genus = "Vibrio" , species = "cholerae" , strain = "MJ" )

all_outgroups = [E_coli,S_enterica,Serratia,P_multocida,P_profundum,V_cholerae]

#clustered_thorsellia = Clustering(all_assemblies, work_root + "MCL_clusters/", "all_thorsellia")
clustered_all = Clustering(all_assemblies + all_outgroups, work_root + "MCL_clusters/", "all_with_outgroups")

clustered_all_nucl = Clustering(all_assemblies + all_outgroups, work_root + "MCL_clusters/all_nucl_w_outgroups/", "all_nucleotide_with_outgroups")

#wolfies = SixteenS(name = "wolbachia_nd_friends", path = work_root + "wolbachia/", raw = work_root + "wolbachia/wolbachia_nd_friends.fasta")
all_wolfies = SixteenS(name = "all_wolfies", path = work_root + "wolbachia/", raw = work_root + "wolbachia/all_wolbachia.fasta")
pap_wolbs = SixteenS("paper_wolbachs",path = work_root + "wolbachia/", raw = work_root + "wolbachia/paper_wolbachias.fasta")

if __name__ == '__main__':
    run()



