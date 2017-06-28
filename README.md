SynerClust README


# Dependencies
- Python-2.7.x
- NetworkX (Python package) http://networkx.github.io/
- Blast+ ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

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
#### Without UGE:
The minimal command to run SynerClust is the following:
<pre><code>path/to/SynerClust/bin/synerclust.py -r path/to/data_catalog.txt -w working/directory/ -t path/to/newick/tree.nwk [-n number_of_cores]</pre></code>

You then need to run the script indicated (all tasks can be run in parallel on the grid):
<pre><code>./genomes/needed_extractions.cmd.txt</pre></code>

You can then start the actual computation (in part parallelizable on the grid):
<pre><code>./jobs.sh</pre></code>

Once all jobs are finished, to have an easy to read output of the clusters, simply run the "post_process_root.sh" script. This will, among others, generate a final_clusters.txt and clusters_to_locus.txt file with the results in the root node.


#### With UGE:
Initialize your environnement (if on UGE):
<pre><code>use Python-2.7
use UGER</pre></code>

The minimal command to run SynerClust is the following:
<pre><code>path/to/SynerClust/bin/synergy.py -r path/to/data_catalog.txt -w working/directory/ -t path/to/newick/tree.nwk [-n number_of_cores]</pre></code>

You then need to run the script indicated:
<pre><code>python path/to/SynerClust/uger_auto_submit_simple.py -f genomes/needed_extractions.cmd.txt -tmp TMP_FOLDER</pre></code>

You can then start the actual computation (parallelizable on the grid):
<pre><code>./uger_jobs.sh</pre></code>

If they are more jobs than your queue allows, run:
<pre><code>path/to/SynerClust/uger_auto_submit.py -f ./uger_jobs.sh -l queue_size_limit</pre></code>

Once all jobs are finished, to have an easy to read output of the clusters, simply run the "post_process_root.sh" script. This will, among others, generate a final_clusters.txt and clusters_to_locus.txt file with the results in the root node.
