#SYNERGY2_PATHWF_ClusterPostProcessing.py #WORKING_DIRgenomes/ #WORKING_DIRnodes/#CHILD1/locus_mappings.pkl
#SYNERGY2_PATHWF_ClusterPostProcessing.py #WORKING_DIRgenomes/ #WORKING_DIRnodes/#CHILD2/locus_mappings.pkl
#SYNERGY2_PATHWF_RunBlast.py #WORKING_DIRnodes/ #NODE #BLAST_EVAL #NUM_CORES #CHILD1 #CHILD2
#SYNERGY2_PATHWF_MakeRoughClusters.py -dir #WORKING_DIRnodes/ -node #NODE -m #MIN_BEST_HIT -F #MIN_SYNTENIC_FRACTION #CHILD1 #CHILD2
#SYNERGY2_PATHWF_RefineClusters.py -dir #WORKING_DIRnodes/ -t #NUM_CORES -node #NODE #CHILD1 #CHILD2 #NOSYNTENY
#SYNERGY2_PATHWF_FinalizeNode.py -dir #WORKING_DIRnodes/#NODE/ -node #NODE -t #NUM_CORES -dist #MUTRATE
