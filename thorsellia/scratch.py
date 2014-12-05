
sari_path = "/home/moritz/people/sari/pelodicton/pelodictyon/"

parvum = Assembly("parvum",  None  , None  , sari_path + "Chlorobaculum_parvum_NCIB_8327.fasta" , genus = "Chlorobaculum" , species = "parvum" , strain = "NCIB_8327" )
limicola = Assembly("limicola",  None  , None  , sari_path + "Chlorobium_limicola_DSM_245.fasta" , genus = "Chlorobium" , species = "limicola" , strain = "DSM_245" )
chromatii = Assembly("chromatii",  None  , None  , sari_path + "Chlorobium_chlorochromatii_CaD3.fasta" , genus = "Chlorobium" , species = "chlorochromatii" , strain = "CaD3" )
luteolum = Assembly("luteolum",  None  , None  , sari_path + "Chlorobium_luteolum_DSM_273.fasta" , genus = "Chlorobium" , species = "luteolum" , strain = "DSM_273" )
BS1 = Assembly("BS1",  None  , None  , sari_path + "Chlorobium_phaeobacteroides_BS1.fasta" , genus = "Chlorobium" , species = "phaeobacteroides" , strain = "BS1" )
DSM_266 = Assembly("DSM_266",  None  , None  , sari_path + "Chlorobium_phaeobacteroides_DSM_266.fasta" , genus = "Chlorobium" , species = "phaeobacteroides" , strain = "DSM_266" )
phaeovibrioides = Assembly("phaeovibrioides",  None  , None  , sari_path + "Chlorobium_phaeovibrioides_DSM_265.fasta" , genus = "Chlorobium" , species = "phaeovibrioides" , strain = "DSM_265" )
tepidum = Assembly("tepidum",  None  , None  , sari_path + "Chlorobium_tepidum_TLS.fasta" , genus = "Chlorobium" , species = "tepidum" , strain = "TLS" )
pelodictyon = Assembly("pelodictyon",  None  , None  , sari_path + "../Contigs.pelodict.fasta" , genus = "Pelodictyon" , species = "phaeoclathratiforme" , strain = "BU-1" )
bin_52_2 = Assembly("bin_52_2",  None  , None  , sari_path + "../../bin_52_2_AKA_TheBin.fasta" , genus = "Pelodictyon" , species = "halsjarvii" , strain = "putative" )
ignavi = Assembly("ignavi",  None  , None  , sari_path + "Ignavibacterium_album_JCM_16511.fasta" , genus = "Ignavivacterium" , species = "album" , strain = "JCM_16511" )

saris = [parvum, limicola, chromatii, luteolum, BS1, DSM_266, phaeovibrioides, tepidum, pelodictyon, bin_52_2]
sari_clusts = Clustering(saris, sari_path + "MCL_clusters/tree_1/", "sari_clusts", seq_type = "genes")


val_path = "/home/moritz/GEFES/views/projects/oilsand/binning/concoct_5000/bins/"
bin_18 = Assembly("bin_18",  None  , None  , val_path + "bin_18/contigs.fasta" , genus = "unknown" , species = "unkown" , strain = "oilsands", type_scaffold = "/home/moritz/people/valerie/OilSands/methano.fasta" )
