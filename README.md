SynerClust README


#Dependencies#
- Python-2.7.x
- NetworkX (Python package) http://networkx.github.io/
- MUSCLE http://www.drive5.com/muscle/downloads.htm
- Blast+ ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

- FastTree http://meta.microbesonline.org/fasttree/#Install  
	For FastTree, it is best to download the C source code and compile it using the following command line :  
<pre><code>gcc -DUSE_DOUBLE -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm </pre></code>


#Installation#
You create a file (commands.install in this example) containing the install settings for ease of use. It should contain the following information in this order :  
<pre><code>/absolute/path/to/FastTree_executable
/absolute/path/to/Muscle_executable
path/to/blast+/bin/   can be left blank if already in the PATH
default e-value threshold
default number of threads to use  </code></pre>
	

You can then install SynerClust by running the following command from the root folder:
<code><pre>python INSTALL.py < commands.install</code></pre>


#Input Data#
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

	
#Running#
####Without UGE:####
The minimal command to run SynerClust is the following:
<pre><code>path/to/SynerClust/bin/WF_runSynergy2.py -r path/to/data_catalog.txt -w working/directory/ -t path/to/newick/tree.nwk</pre></code>

You can optionally add settings as in this example with default values:
<pre><code>-a 0.01 -b 1.0 -g 5.0 -G 0.05 -L 0.05 -D 1.2</pre></code>

You then need to run the script indicated (all tasks can be run in parallel on the grid):
<pre><code>./genomes/needed_extractions.cmd.txt</pre></code>

Then re-run the first command
<pre><code>path/to/SynerClust/bin/WF_runSynergy2.py -r path/to/data_catalog.txt -w working/directory/ -t path/to/newick/tree.nwk</pre></code>

You can then start the actual computation (in part parallelizable on the grid):
<pre><code>./jobs.sh</pre></code>

Once all jobs are finished, to have an easy to read output of the clusters, run:
<pre><code>path/to/SynerClust/bin/ClusterPostProcessing.py genomes/ nodes/N_**********_***************/locus_mappings.pkl n</pre></code>
Where the * are the id of the node (for the root, the first number is the depth of the tree, which is also the highest number generated on this run) and n the number of genomes included at that node.
This will generate a final_clusters.txt and clusters_to_locus.txt file with the results.


####With UGE:####
Initialize your environnement (if on UGE):
<pre><code>use Python-2.7
use UGER</pre></code>

The minimal command to run SynerClust is the following:
<pre><code>path/to/SynerClust/bin/WF_runSynergy2.py -r path/to/data_catalog.txt -w working/directory/ -t path/to/newick/tree.nwk</pre></code>

You can optionally add settings as in this example with default values:
<pre><code>-a 0.01 -b 1.0 -g 5.0 -G 0.05 -L 0.05 -D 1.2</pre></code>

You then need to run the script indicated:
<pre><code>python path/to/SynerClust/uger_auto_submit_simple.py -f genomes/needed_extractions.cmd.txt -tmp TMP_FOLDER</pre></code>

Then re-run the first command
<pre><code>path/to/SynerClust/bin/WF_runSynergy2.py -r path/to/data_catalog.txt -w working/directory/ -t path/to/newick/tree.nwk</pre></code>

You can then start the actual computation (in part parallelizable on the grid):
<pre><code>./uger_jobs.sh</pre></code>

If they are more jobs than your queue allows, run:
<pre><code>path/to/SynerClust/uger_auto_submit.py -f ./uger_jobs.sh -l queue_size_limit</pre></code>

Once all jobs are finished, to have an easy to read output of the clusters, run:
<pre><code>path/to/SynerClust/bin/ClusterPostProcessing.py genomes/ nodes/N_**********_***************/locus_mappings.pkl n</pre></code>
Where the * are the id of the node (for the root, the first number is the depth of the tree, which is also the highest number generated on this run) and n the number of genomes included at that node.
This will generate a final_clusters.txt and clusters_to_locus.txt file with the results.
