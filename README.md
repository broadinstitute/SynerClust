SynerClust README


# Dependencies
- Python-2.7.x
- NumPy (Python package) https://www.scipy.org/scipylib/download.html
- NetworkX (Python package) https://pypi.python.org/pypi/networkx/
- Blast+ (tested with v2.6.0) ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

Already included:
- MUSCLE http://www.drive5.com/muscle/downloads.htm
- FastTree http://meta.microbesonline.org/fasttree/#Install  


# Installation
You can install SynerClust by running the following command from the main folder:
<code><pre>python INSTALL.py</code></pre>

The default considers that Blast+ is in your path. If that is not the case, please use the "-e" option to specify the path to the Blast+ bin folder.


# Input Data
data_catalog.txt should be formatted as the following example (paths can be relative of absolute paths):
<code><pre>
//
Genome	Esch_coli_H296
Sequence	Esch_coli_H296/Esch_coli_H296.genome
Annotation	Esch_coli_H296/Esch_coli_H296_PRODIGAL_2.annotation.gff3
//
Genome	Esch_coli_H378_V1
Sequence	Esch_coli_H378_V1/Esch_coli_H378_V1.genome
Annotation	Esch_coli_H378_V1/Esch_coli_H378_V1_PRODIGAL_2.annotation.gff3
//
</code></pre>

	
# Running
### On a single machine:
The minimal command to run SynerClust is the following:
<pre><code>/path/to/SynerClust/bin/synerclust.py -r /path/to/data_catalog.txt -w /working/directory/ -t /path/to/newick/tree.nwk [-n number_of_cores] [--run single]</pre></code>

If you use the option "--run single" that is all you need to do! The Results will be written in the root folder but symlinks can be found in in the main "results" folder.


If you prefer to run step by step, the next steps are:

Run the script indicated (all tasks can be run in parallel on a grid):
<pre><code>/working/directory/genomes/needed_extractions.cmd.txt</pre></code>

You can then start the actual computation (in part parallelizable on a grid):
<pre><code>/working/directory/jobs.sh</pre></code>

Once all jobs are finished, to have an easy to read output of the clusters, simply run:
<pre><code>/working/directory/post_process_root.sh</pre></code>
This will, among others, generate a final_clusters.txt and clusters_to_locus.txt file with the results in the "results" folder (linking to the root node).


### With an UGE cluster:
Initialize your environnement (if on UGE):
<pre><code>use Python-2.7
use UGER</pre></code>

The minimal command to run SynerClust is the following:
<pre><code>/path/to/SynerClust/bin/synerclust.py -r /path/to/data_catalog.txt -w /working/directory/ -t /path/to/newick/tree.nwk [-n number_of_cores] [--run uger]</pre></code>

If you use the option "--run uger" that is all you need to do! The Results will be written in the root folder but symlinks can be found in in the main "results" folder.


If you prefer to run step by step, the next steps are:

Run the script indicated:
<pre><code>python /path/to/SynerClust/uger_auto_submit_simple.py -f /working/directory/genomes/needed_extractions.cmd.txt -tmp TMP_FOLDER</pre></code>

You can then start the actual computation (parallelizable on the grid):
<pre><code>/working/directory/uger_jobs.sh</pre></code>

If they are more jobs than your queue allows, run:
<pre><code>/path/to/SynerClust/uger_auto_submit.py -f /working/directory/uger_jobs.sh -l queue_size_limit [-n number_of_cores_per_job]</pre></code>

Once all jobs are finished, to have an easy to read output of the clusters, simply run the 
<pre><code>/working/directory/post_process_root.sh</pre></code>
This will, among others, generate a final_clusters.txt and clusters_to_locus.txt file with the results in the "results" folder (linking to the root node).


# Help/Questions

### Output files
clusters.txt: Contains a list of all genes grouped by clusters with annotation information. One gene per line, clusters are seperated by empty lines.

final_clusters.txt: List of all the transcript IDs per cluster and count of genes and taxas. One cluster per line.

clust_to_trans.txt: List of all transcripts IDs, their cluster assignation and the most commonly encountered gene name that is not "hypothetical protein". One gene per line.

cluster_dist_per_genome.txt: Table of number of genes from each taxas per cluster. One cluster per line.

### Keeping track of progress

A file "completion.txt" is generated in the main folder of the run where you can keep track of what nodes have finished being processed.

In case a node encountered an error, another file, "not_completed.txt" is generated containing the node identifier.

### Running SynerClust on an extended dataset

Simply rerun SynerClust as you did on the firsrt dataset, but set the working directory to the same one as the first run.

### List of Parameters and their meaning
-t SPECIES_TREE, --tree SPECIES_TREE  
&nbsp;&nbsp;&nbsp;&nbsp;Species tree relating all of the genomes to be analyzed. (Required)

-r COBRA_REPO, --repo COBRA_REPO  
&nbsp;&nbsp;&nbsp;&nbsp;Complete path to data_catalog in the repository containing your genomic data. (Required)

-w WORKING_DIR, --working WORKING_DIR  
&nbsp;&nbsp;&nbsp;&nbsp;Complete path to the working directory for this analysis. (Required)

-m MIN_BEST_HIT, --min_best_hit MIN_BEST_HIT  
&nbsp;&nbsp;&nbsp;&nbsp;Minimal % of match length for Blastp hits compared to best one. (default = 0.8)

-B BLAST_EVAL, --blast_eval BLAST_EVAL  
&nbsp;&nbsp;&nbsp;&nbsp;Minimal e-value for Blastp hits. (default = 0.00001)

-l LOCUS_FILE, --locus LOCUS_FILE  
&nbsp;&nbsp;&nbsp;&nbsp;A locus_tag_file.txt that corresponds to the data in this repository

-N CODED_NWK_FILE, --newick_tag CODED_NWK_FILE  
&nbsp;&nbsp;&nbsp;&nbsp;Output file for the newick tree using tag names and number of genomes as distances.

-n NUM_CORES, --num_cores NUM_CORES  
&nbsp;&nbsp;&nbsp;&nbsp;The number of cores used for blast analysis (-a flag), (default = 4)

-F MINSYNFRAC, --min_syntenic_fraction MINSYNFRAC  
&nbsp;&nbsp;&nbsp;&nbsp;Minimum common syntenic fraction required for two genes from the same species to be considered paralogs, range [0.0,1.0], default=0.7

-D DIST, --dist DIST  
&nbsp;&nbsp;&nbsp;&nbsp;Maximum FastTree distance between a representative sequence and sequences being represented for representative selection. (default = 1.2)

-s SYNTENY_WINDOW, --synteny_window SYNTENY_WINDOW  
&nbsp;&nbsp;&nbsp;&nbsp;Distance in base pairs that will contribute to upstream and downstream to syntenic fraction. The total window size is [int]*2. (default = 6000)

--no-synteny  
&nbsp;&nbsp;&nbsp;&nbsp;Disable use of synteny (required if information not available).

  --run {none,single,uger}  
&nbsp;&nbsp;&nbsp;&nbsp;Specify if you want all computation to be run directly. Use "single" to run the local machine or "uger" to submit to the UGE grid.

  --alignment {none,scc,all}  
&nbsp;&nbsp;&nbsp;&nbsp;Specify if you want cluster alignments using MUSCLE to be computed and written for the root node. Use "all" if you want all clusters to be aligned or "scc" if you only want Single Copy Core clusters to be aligned.


