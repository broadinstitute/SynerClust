#! /bin/bash

#$ -cwd
#$ -q long
#$ -P gscid

#$ -l m_mem_free=2g
#$ -e /cil/shed/sandboxes/cgeorges/uger/error.err
#$ -o /cil/shed/sandboxes/cgeorges/uger/out.log

source /broad/software/scripts/useuse
reuse -q Python-2.7 BLAST+

#SYNERGY2_PATHWF_ClusterPostProcessing.py #WORKING_DIRgenomes/ #WORKING_DIRnodes/#CHILD1/locus_mappings.pkl
#SYNERGY2_PATHWF_ClusterPostProcessing.py #WORKING_DIRgenomes/ #WORKING_DIRnodes/#CHILD2/locus_mappings.pkl
#SYNERGY2_PATHWF_RunBlast.py #WORKING_DIRnodes/ #NODE #BLAST_EVAL #NUM_CORES #CHILD1 #CHILD2
#SYNERGY2_PATHWF_MakeRoughClusters.py -dir #WORKING_DIRnodes/ -node #NODE -m #MIN_BEST_HIT -F #MIN_SYNTENIC_FRACTION #CHILD1 #CHILD2
#SYNERGY2_PATHWF_RefineClusters_leaf_centroid_newmatrix.py -dir #WORKING_DIRnodes/ -node #NODE -alpha #ALPHA -beta #BETA -gamma #GAMMA -gain #GAIN -loss #LOSS #CHILD1 #CHILD2
#SYNERGY2_PATHWF_FinalizeNode_threaded.py -dir #WORKING_DIRnodes/#NODE/ -node #NODE -t #NUM_CORES -dist #HAMMING
