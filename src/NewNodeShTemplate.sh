#! /bin/bash

#$ -cwd
#$ -q long
#$ -P gscid

#$ -l m_mem_free=2g
#$ -e /cil/shed/sandboxes/cgeorges/uger/error.err
#$ -o /cil/shed/sandboxes/cgeorges/uger/out.log

source /broad/software/scripts/useuse
reuse -q Python-2.7 Java-1.7

#SYNERGY2_PATHWF_ClusterPostProcessing.py #WORKING_DIRgenomes/ #WORKING_DIRnodes/#CHILD1/locus_mappings.pkl
#SYNERGY2_PATHWF_ClusterPostProcessing.py #WORKING_DIRgenomes/ #WORKING_DIRnodes/#CHILD2/locus_mappings.pkl
#SYNERGY2_PATHWF_RunBlast.py #WORKING_DIRnodes/ #NODE #BLAST_EVAL #NUM_CORES #CHILD1 #CHILD2
#SYNERGY2_PATHWF_MakeRoughClusters.py #WORKING_DIRnodes/ #NODE #MIN_BEST_HIT #HOMOLOGY_SCALE #SYNTENY_SCALE #NUM_HITS #MIN_SYNTENIC_FRACTION #CHILD1 #CHILD2
#SYNERGY2_PATHWF_RefineClusters_leaf_centroid_newmatrix.py #WORKING_DIRnodes/ #ID.5 #CMDS_PER_JOB #NODE #ALPHA #GAMMA #GAIN #LOSS #CHILD1 #CHILD2
#SYNERGY2_PATHWF_FinalizeNode_threaded.py #WORKING_DIRnodes/#NODE/ #NODE #NUM_CORES #HAMMING
