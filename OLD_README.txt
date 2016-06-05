README

OVERVIEW
Synergy2 has 3 major phases.  First, a file system is established from gff3 and fasta files, which can be used for every subsequent Synergy2 analysis performed on the same set of genomes. A species tree is either provided by the user or generated using our pipeline based on AMPHORA (see 3.2.1)


1 Dependencies
1.1 Files
1.1.1 Genome-specific files (GFF3 structural annotation and fasta-formatted assembly)
1.1.2 Species tree

1.2 Programs
1.2.1 QuickTree
1.2.2 Muscle
1.2.3 Blast

1.3 Packages
1.3.1 WorkFlow
A functional instance of WorkFlow (available at http://sourceforge.net/projects/tigr-workflow/) is required to run Synergy2.
You can install WorkFlow to run using a grid, or to run linearly.  Synergy2 only requires linear execution.
1.3.2 Python 2.6
1.3.2.1 Standard: sys, string, random, os, numpy, scipy, pickle, re, math, getopt, subprocess, time
1.3.2.2 Additional: NetworkX (http://networkx.lanl.gov/index.html)

2 Installation
To install Synergy2, extract the tarball in the permanent location of the program.  From within that directory, run INSTALL.py.
2.1 Paths to be set
QuickTree
Muscle
Blast
2.2 Other Variables - these are your chosen defaults, and may be changed at run-time of a given Synergy2 instance
Number of cores to be used by BLAST
E-value used by BLAST



USAGE

3 File System Set-up (satisfying file dependency requirements)
3.1 Genomic data extraction
3.1.1 Data Specifications File
A file containing the complete path to the genome sequence (fasta format) and corresponding structural genome annotation (GFF3 format, http://www.sequenceontology.org/gff3.shtml) must be provided.  The sequence ID (1st column) of the GFF3 file must match exactly to a sequence header in the fasta file. The 9th column of the GFF3 must contain "ID=<string>", where <string> is a unique identifier of that particular gene.  See <EXAMPLE> for this file format.
//
Genome Bsub_E342_V1
Sequence /path/to/bsub_sequence.fa
Assembly /path/to/bsub_annot.gff3
//

3.2 Species tree
The species tree must be unrooted and in Newick format.  Branch lengths may be in decimal or scientific notation. Bootstrap values, while allowed, will not factor into the algorithm.  Branch lengths will only be used to root the tree by its midpoint. If the tree has ancestral nodes with more than 2 children, the children will be randomly grouped into pairs. 

If the resolution of the initial species tree is not satisfactory, a new species tree can be generated from the single copy core after Synergy2 has been run.  At that point, Synergy2 may be re-run with the updated tree.  Computations will only be performed where the topology has changed.

3.2.1 Using AMPHORA2
If no species tree is available, the user can create a species tree based on 39 protein HMMs defined by AMPHORA2.  (http://wolbachia.biology.virginia.edu/WuLab/Software.html)


4 Initializing Synergy2
After all of the data requirements have been satisfied, the user runs WF_runSynergy2.py.  This script defines the WorkFlow instance and creates all of the dependent files.  This is also the control point for the user to define the variable parameters.
4.1 Variables
4.1.1 Computational Resources
BLAST - # of cores, default=4:
This is configured during install, but can be modified here.

Wait for File - update frequency, default=60s:
How often, in seconds, a sub workflow checks to see if its dependencies have been satisfied. Generally should not be changed.

4.1.2 Results
These variables impact the output of Synergy2. 
4.1.2.1 Rough Clusters

Blast e-value, default=0.01:
Determines the minimum quality of acceptable hit.  Only best hits will be considered to define rough clusters, but all hits will be used to compute syntenic fractions and distance matrices. Consequently, this cut-off is not stringent.

Syntenic window size, default=5kb:
Extends to the left and right of the gene in question.  A 5kb window evaluates a total of 10kb surrounding each gene. Raise or lower based on the assumed range of syntenic conservation in your analysis.

Minimum "Best Hit" score, default=0.5: 
This determines the tightness of the rough clusters.  A lower value creates larger rough clusters as it allows for more variation, while a higher value does the opposite.  This score is representative of the quality and coverage of the hit. The coverage is calculated with respect the percent identity and query vs target coverage of the hit.

Scaling of contribution for graph edges, default homology=0.01, default synteny=1.0:
Impacts the weight of homology (BLAST) and synteny when a cluster's gene tree is rooted or evaluated. Because homology is used to create the partitions for the gene tree, it is scaled back.  If the genomes are particularly fragmented, homology weight may be increased to allow for differentiation.

4.1.2.2 Tree Rooting

gamma, default=10.0
A factor in the gene tree rooting equation that affects the impact of gain and loss events

gain, default=0.05 
The probability of a gain event at a given ancestral node 

loss, default=0.05 
The probability of a loss event at a given ancestral node 

5 Running the WorkFlow instance 
WF_runSynergy2.py returns a list of commands that need to be executed in order to set up the working environment.  These commands, found in genomes/needed_extractions.cmd.txt, can be run serially or distributed to a grid. After the working environment has been set up, re-run WF_runSynergy2.py the same as before. This will create a [flow].xml and a [flow].ini file in your Synergy2 working directory, where [flow] is defined by -f or --flow.  You are now ready to run the actual algorithm.  

Set up the command line environment.

source /path/to/workflow/exec_env.tcsh

use LSF (make sure your grid distribution module is loaded)

Launch instance_commands.txt to the grid.  If you have queues for different time periods, an average Synergy2 run of 150 E.coli genomes can be completed in about 12 hours across 50 nodes.

Then, the WorkFlow instance may be dispatched as follows:

RunWorkflow -t [flow].xml -c [flow].ini -i [instance].xml &

This is not necessary for the execution of Synergy2, but will allow you to monitor its progress.