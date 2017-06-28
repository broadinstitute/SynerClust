cat genomes/locus_tag_file.txt | grep N_0 | awk '{print "#SYNERCLUST_PATHClusterPostProcessing.py genomes/ nodes/"$2"/locus_mappings.pkl "$3}' > post_process_all.sh
chmod +x post_process_all.sh

