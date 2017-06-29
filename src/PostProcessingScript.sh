#!/usr/bin/env bash

echo "#!/usr/bin/env bash" > post_process_all.sh
echo "" >> post_process_all.sh
cat genomes/locus_tag_file.txt | grep -vP "\t1$" | grep N_0 | awk '{print "#SYNERCLUST_PATHClusterPostProcessing.py genomes/ nodes/"$2"/locus_mappings.pkl "$3" -a '$1'"}' >> post_process_all.sh
chmod +x post_process_all.sh

echo "#!/usr/bin/env bash" > post_process_root.sh
echo "" >> post_process_root.sh
sort -nk4,4 post_process_all.sh | tail -1 >> post_process_root.sh
chmod +x post_process_root.sh
