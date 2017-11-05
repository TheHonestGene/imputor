import pytest
from imputor.core import genotype_parser

class TestPrepare:

    def test_convert_csv_to_hdf5(self, genotype_as_csv, genotype_as_hdf5, tmpdir):
        # TODO call genotype_parser and test that output hdf5 is correct (compare with genotype_as_hdf5)
        raise Exception('Test not implemented')

    def test_convert_genotype(self, genotype_as_hdf5, nt_map, genotype_nt_as_hdf5, tmpdir ):
        # TODO call genotype_parser and test that output hdf5 is correct (compare with genotype_nt_as_hdf5)
        raise Exception('Test not implemented')