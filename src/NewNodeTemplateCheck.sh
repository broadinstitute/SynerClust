#!/usr/bin/env bash

#SYNERCLUST_PATHWF_ClusterPostProcessing.py #WORKING_DIRgenomes/ #WORKING_DIRnodes/#CHILD1/locus_mappings.pkl
#SYNERCLUST_PATHWF_ClusterPostProcessing.py #WORKING_DIRgenomes/ #WORKING_DIRnodes/#CHILD2/locus_mappings.pkl
#SYNERCLUST_PATHWF_RunBlast.py #WORKING_DIRnodes/ #NODE #BLAST_EVAL #NUM_CORES #CHILD1 #CHILD2
#SYNERCLUST_PATHWF_MakeRoughClusters.py -dir #WORKING_DIRnodes/ -node #NODE -m #MIN_BEST_HIT -F #MIN_SYNTENIC_FRACTION -i #MIN_PERCENT_IDENTITY -c #MIN_MATCH_COVERAGE #CHILD1 #CHILD2
#SYNERCLUST_PATHWF_RefineClusters.py -dir #WORKING_DIRnodes/ -t #NUM_CORES -node #NODE #CHILD1 #CHILD2 #NOSYNTENY
#SYNERCLUST_PATHWF_FinalizeNode.py -dir #WORKING_DIRnodes/#NODE/ -node #NODE -t #NUM_CORES -dist #MUTRATE
if [ -f #WORKING_DIRnodes/#NODE/NODE_COMPLETE ]; then echo "#NODE completed" >> #WORKING_DIRcompletion.txt ; else echo "#NODE did not complete" >> not_completed ; fi
