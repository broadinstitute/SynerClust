#! /bin/bash

#$ -cwd
#$ -q short
#$ -P gscid

#$ -l h_vmem=2g
#$ -binding linear:#NUM_CORES
#$ -pe smp #NUM_CORES
#$ -e #WORKING_DIRerror.err
#$ -o #WORKING_DIRout.log

source /broad/software/scripts/useuse
reuse -q Python-2.7 BLAST+ Coreutils

#SYNERGY2_PATHWF_ClusterPostProcessing.py #WORKING_DIRgenomes/ #WORKING_DIRnodes/#CHILD1/locus_mappings.pkl
#SYNERGY2_PATHWF_ClusterPostProcessing.py #WORKING_DIRgenomes/ #WORKING_DIRnodes/#CHILD2/locus_mappings.pkl
#SYNERGY2_PATHWF_RunBlast.py #WORKING_DIRnodes/ #NODE #BLAST_EVAL #NUM_CORES #CHILD1 #CHILD2
#SYNERGY2_PATHWF_MakeRoughClusters.py -dir #WORKING_DIRnodes/ -node #NODE -m #MIN_BEST_HIT -F #MIN_SYNTENIC_FRACTION #CHILD1 #CHILD2
#SYNERGY2_PATHWF_RefineClusters_leaf_centroid_newmatrix_bis_threaded.py -dir #WORKING_DIRnodes/ -t #NUM_CORES -node #NODE #CHILD1 #CHILD2 #NOSYNTENY
#SYNERGY2_PATHWF_FinalizeNode_threaded.py -dir #WORKING_DIRnodes/#NODE/ -node #NODE -t #NUM_CORES -dist #MUTRATE
