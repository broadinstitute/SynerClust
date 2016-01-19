EXAMPLE DIRECTIONS

1. cd /Synergy2/example/
2. edit data_catalog.txt to reflect the complete path to each file
~ $ pwd
/home/utilities/Synergy2/example/
//
Genome	Esch_coli_H296
Sequence	/home/utilities/Synergy2/example/Esch_coli_H296/Esch_coli_H296.genome
Annotation	/home/utilities/Synergy2/example/Esch_coli_H296/Esch_coli_H296_PRODIGAL_2.annotation.gff3
//
Genome	Esch_coli_H378_V1
Sequence	/home/utilities/Synergy2/example/Esch_coli_H378_V1/Esch_coli_H378_V1.genome
Annotation	/home/utilities/Synergy2/example/Esch_coli_H378_V1/Esch_coli_H378_V1_PRODIGAL_2.annotation.gff3
//

3. A species tree has been provided (tree.nwk). It is unrooted and in newick format. If a species tree did not exist, one could be constructed with AMPHORA2.
4. Create a working directory for this clustering analysis
example $ mkdir test
5. Move to the test directory
example $ cd test/
6. At this point, all data dependencies have been satisfied and the Synergy2 working directory can be populated with data. For this example, the default values are used.  To see the current defaults:
test $ ../../bin/WF_runSynergy2.py --help

7. Initialize data structures
test $../../bin/WF_runSynergy2.py -t ../tree.nwk -r ../data_catalog.txt -w /home/utilities/Synergy2/example/test/ -f test_run
Wrote locus tags to locus_tag_file.txt
Launch this command file on the grid: /home/utilities/Synergy2/example/test/genomes/needed_extractions.cmd.txt

8. Run commands in needed_extractions.cmd.txt

9. Re-run the command from step 7
test $../../bin/WF_runSynergy2.py -t ../tree.nwk -r ../data_catalog.txt -w /home/utilities/Synergy2/example/test/ -f test_run
Wrote locus tags to locus_tag_file.txt
init tree lib
reading genome to locus
reading tree
((Esch_coli_H378_V1:0.001,Esch_coli_TA014_V1:0.0014):0.002,((Esch_coli_H461_V1:0.0012,Esch_coli_R527_V1:0.001):0.0015,Esch_coli_H296:0.003):0.0025);
parsing tree
VDY;ZGN XSL
IBX;OXY IGQ
IGQ;SVM HAY
XSL:0.002,HAY:0.0025
rooting
('XSL', 'HAY')
XSL 2 [('XSL', 'VDY'), ('XSL', 'ZGN')]
edge from root to XSL
XSL 3 {'VDY': {'weight': 0.0014}, 'root': {'weight': 0.001}, 'ZGN': {'weight': 0.001}}
HAY 2 [('HAY', 'SVM'), ('HAY', 'IGQ')]
edge from root to HAY
HAY 3 {'SVM': {'weight': 0.0030000000000000001}, 'root': {'weight': 0.001}, 'IGQ': {'weight': 0.0015}}
HAY;XSL RPJ
initializing
dependicizing
set number  1 IGQ
set number  2 XSL
{1: [('1.1.1', 'IGQ'), ('1.1.2', 'HAY'), ('1.1.3', 'root')], 2: [('1.2.1', 'XSL')]}
Wrote locus tags to locus_tag_file.txt

NOTE: The 3 letter strings following each genome will not be the same for your example.  They are stored in 'genomes/locus_tag_file.txt', and therefore will be consistent between all runs executed in the CWD.

10. Make sure environment is properly initialized
test $ use GridEngine (or use LSF or whatever grid you are using)
test $ source /path/to/your/workflow/install/exec_env.tcsh

11. Launch instance_commands.txt to the grid.  This is where the bulk of the computation occurs.

12. To monitor the progress of the workflow, launch a Workflow instance:
test $ RunWorkflow -t test_run.xml -c test_run.ini -i test_instance.xml &

To monitor progress in a GUI environment, make sure you are forwarding X11 packets, then:
test $ MonitorWorkflow -i test_instance.xml &

To monitor progress in your terminal:
test $ ../../WF_monitor.py test_instance.xml 30
NOTE: 30 is the number of seconds between refresh


13. After WorkFlow instance has completed, the results need to be translated back to identifiers that correspond to the user-specified data.  This can be done for any node, but to get the clustering data for all genomes, it must be done for the root node, which is always test/nodes/root/. Run ClusterPostProcessing.py.
test $ ../../bin/ClusterPostProcessing.py genomes/ nodes/root/locus_mapping.pkl 5
total genes:  24570
pairs: 43456
scc: 3854
mcc: 51
aux: 1280
orphans: 1270
non-orphan clusters: 5185

This will generate several files in nodes/root/.
final_clusters.txt - each line of this file describes 1 cluster
clust_to_trans.txt - a 2-column mapping of cluster ID to gene ID
tuple_pairs.pkl - used to calculate the Jaccard Coefficient between two clustering analyses with JaccardPickleComparator.py
summary_stats.txt - used as input for SummaryStatsParse.py

To see if the results you have generated match the results in the example:
test $ ../../bin/JaccardPickleComparator.py ../results/tuple_pairs.pkl nodes/root/tuple_pairs.pkl
OR
test $ ../../SummaryStatsParse.py nodes/root/summary_stats.txt 5
