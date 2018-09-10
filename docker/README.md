# Run SynerClust via Docker

	% docker run --rm -v `pwd`:`pwd` synerclust/synerclust synerclust.py \
	-r `pwd`/path/to/data_catalog.txt \
	-w `pwd`/working/directory/ \
	-t `pwd`/path/to/newick/tree.nwk \
	-n 4	
	
