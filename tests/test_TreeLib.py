import TreeLib
from nose.tools import assert_equals
from mock import MagicMock
from mock import patch


class Test_TreeLib:
	def setup(self):
		mock_open = MagicMock()
		with patch('__builtin__.open', mock_open):
			manager = mock_open.return_value.__enter__.return_value
			manager.read.return_value = "(A:0.1,B:0.2,(C:0.3,D:0.4)5.51:0.5);"

			self.tree = TreeLib.Tree("./example/tree.nwk", "./tests/data/genomes/locus_tag_file.txt")
			self.tree.readTree()

	def teardown(self):
		self.tree = None

	def test_readTree(self):
		assert_equals(self.tree.tree_string, "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);")

