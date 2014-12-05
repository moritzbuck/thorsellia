import sh
import os
from Bio import SeqIO
from thorsellia.parameters import *
from pandas import DataFrame

blast_head = ["query","subject","identity","length","mismatches", "gapopenings", "querystart","queryend","subjectstart","subjectend","Evalue","bitscore"]

blast_dbs = {
    "COG0012": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0012.fna", "sets": ["SCG", "cog"]},
    "COG0016": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0016.fna", "sets": ["SCG", "cog"]},
    "COG0018": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0018.fna", "sets": ["SCG", "cog"]},
    "COG0048": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0048.fna", "sets": ["SCG", "cog"]},
    "COG0049": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0049.fna", "sets": ["SCG", "cog"]},
    "COG0052": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0052.fna", "sets": ["SCG", "cog"]},
    "COG0080": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0080.fna", "sets": ["SCG", "cog"]},
    "COG0081": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0081.fna", "sets": ["SCG", "cog"]},
    "COG0085": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0085.fna", "sets": ["SCG", "cog"]},
    "COG0087": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0087.fna", "sets": ["SCG", "cog"]},
    "COG0088": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0088.fna", "sets": ["SCG", "cog"]},
    "COG0090": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0090.fna", "sets": ["SCG", "cog"]},
    "COG0091": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0091.fna", "sets": ["SCG", "cog"]},
    "COG0092": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0092.fna", "sets": ["SCG", "cog"]},
    "COG0093": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0093.fna", "sets": ["SCG", "cog"]},
    "COG0094": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0094.fna", "sets": ["SCG", "cog"]},
    "COG0096": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0096.fna", "sets": ["SCG", "cog"]},
    "COG0097": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0097.fna", "sets": ["SCG", "cog"]},
    "COG0098": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0098.fna", "sets": ["SCG", "cog"]},
    "COG0099": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0099.fna", "sets": ["SCG", "cog"]},
    "COG0100": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0100.fna", "sets": ["SCG", "cog"]},
    "COG0102": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0102.fna", "sets": ["SCG", "cog"]},
    "COG0103": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0103.fna", "sets": ["SCG", "cog"]},
    "COG0124": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0124.fna", "sets": ["SCG", "cog"]},
    "COG0172": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0172.fna", "sets": ["SCG", "cog"]},
    "COG0184": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0184.fna", "sets": ["SCG", "cog"]},
    "COG0185": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0185.fna", "sets": ["SCG", "cog"]},
    "COG0186": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0186.fna", "sets": ["SCG", "cog"]},
    "COG0197": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0197.fna", "sets": ["SCG", "cog"]},
    "COG0200": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0200.fna", "sets": ["SCG", "cog"]},
    "COG0201": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0201.fna", "sets": ["SCG", "cog"]},
    "COG0202": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0202.fna", "sets": ["SCG", "cog"]},
    "COG0215": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0215.fna", "sets": ["SCG", "cog"]},
    "COG0256": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0256.fna", "sets": ["SCG", "cog"]},
    "COG0495": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0495.fna", "sets": ["SCG", "cog"]},
    "COG0522": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0522.fna", "sets": ["SCG", "cog"]},
    "COG0525": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0525.fna", "sets": ["SCG", "cog"]},
    "COG0533": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0533.fna", "sets": ["SCG", "cog"]},
    "COG0541": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0541.fna", "sets": ["SCG", "cog"]},
    "COG0552": {"path": "/home/moritz/glob/data/blast.dbs/40scg_from_specI/COG0552.fna", "sets": ["SCG", "cog"]}
    }



def make_db(assemblies ,path , seq_type = "genes" ):
    #seq_type can be "genes", "genomes", or "proteins"
    print "Making a", seq_type, "blast db"
    seqs = []
    for a in assemblies:

        fila = {"genes" : a.genes , "genomes" : a.clean_built , "proteins" : a.proteins}[seq_type]

        with open(fila) as assembly:
            tseqs = SeqIO.parse(assembly,"fasta")
            for s in tseqs:
                if seq_type == "genomes" : 
                    s.id = a.name + "_" + s.id
                seqs.append(s)
    with open(path,"w") as combined_fasta:
         SeqIO.write(seqs, combined_fasta, "fasta")
         
    sh.makeblastdb("-in", path, "-dbtype" , "prot" if seq_type == "proteins" else "nucl")

########################################################################################################3
    
def blast(query,output_path, db, alg = "blastp", clean=True,post_process = True, outfmt = 6, evalue = 0.001, word_size = None):
    itype = "gene" if alg == "blastn" or alg == "tblastx" else "protein"
    otype = "protein" if alg == "blastp" or alg == "tblastn" else "gene"
    print "Blasting", query, " to the", otype, "db"
    blasty = sh.Command(alg)


    if not word_size:
        if itype == "gene": word_size = 11
        else : word_size = 3

    blasty("-evalue", evalue, "-db", db, "-query",query, "-outfmt", outfmt , "-word_size", word_size, "-num_threads", threads, _out =  output_path + ".raw.blast")

    if(os.stat(output_path+".raw.blast").st_size > 0 and post_process):
        raw_blast = DataFrame.from_csv(output_path+".raw.blast",header=-1,index_col=[0,1],sep="\t")
        raw_blast.columns = blast_head[2:12]
        raw_blast.index.names = blast_head[0:2]
        asse = ["_".join(b[1].split("_")[:-1]) for b in raw_blast.index]
        ref = [b[1] for b in raw_blast.index]
        query = [b[0] for b in raw_blast.index]
        eids = ",".join([i.split("|")[1] for i in query])
        handle = Entrez.efetch(db= itype, id = eids, retmode="text", rettype="gb")
        descr = ["|".join([r.description,"length:" + str(len(r.seq))]) for r in SeqIO.parse(handle,"genbank")]
        query = ["".join(zizi) for zizi in zip(query,descr)]
        index = MultiIndex.from_arrays([query,asse],names = ["query","assembly"])
    
        data = append(matrix(ref).transpose(),matrix(raw_blast[["identity","length","Evalue","bitscore"]]),1)
        raw_blast = DataFrame(data,index=index, columns = ["gene_id","identity","length","Evalue","bitscore"] )
        processed_blast = raw_blast.groupby(level = ["query","assembly"]).apply(lambda x: x.iloc[ x['length'].argmax() ])
    else:
        raw_blast = DataFrame.from_csv(output_path+".raw.blast",header=-1,index_col=[0,1],sep="\t")
        raw_blast.columns = blast_head[2:12]
        raw_blast.index.names = blast_head[0:2]
        processed_blast = raw_blast
        
    if clean: os.remove(output_path+".raw.blast")
    return processed_blast

                
########################################################################################################3
