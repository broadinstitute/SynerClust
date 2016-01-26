import COBRA_Repo_Handling
from nose.tools import assert_equals
from mock import MagicMock
from mock import patch

class Test_COBRA_Repo_Handling:
	# def setup(self):
	# 	self.rep = RepoParse("../example/data_catalog.txt")

	def setup(self):
		mock_open = MagicMock()
		with patch('__builtin__.open', mock_open):
			manager = mock_open.return_value.__enter__.return_value
			manager.readline.return_value = ["//\n", "Genomes\tEsch_coli_H296\n", "Sequences\t./example/Esch_coli_H296/Esch_coli_H296.genome\n", "Annotations\t./example/Esch_coli_H296/Esch_coli_H296_PRODIGAL_2.annotation.gff3\n", "//\n"]

			self.rep = COBRA_Repo_Handling.RepoParse("../example/data_catalog.txt")
			self.rep.parseRepoFile("../example/")


	def teardown(self):
		print "teardown"

	def test_parseRepoFile(self):
		expected = COBRA_Repo_Handling.Genome("Esch_coli_H296", "./example/Esch_coli_H296/Esch_coli_H296_PRODIGAL_2.annotation.gff3", "./example/Esch_coli_H296/Esch_coli_H296.genome", "null")
		print self.rep.genomes
		assert_equals(self.rep.genomes[0].genome, expected.genome)
		assert_equals(self.rep.genomes[0].locus, expected.locus)
		assert_equals(self.rep.genomes[0].annotation, expected.annotation)
		assert_equals(self.rep.genomes[0].sequence, expected.sequence)
		assert_equals(self.rep.genomes[0].peptide, expected.peptide)
		assert_equals(self.rep.genomes[0].directory, expected.directory)
