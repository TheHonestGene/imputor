from os import path
import pytest
import json


data_path = path.join(path.dirname(__file__), 'data')

@pytest.fixture
def genotype_as_hdf5():
    return path.join(data_path, 'genotype.hdf5')

@pytest.fixture
def genotype_nt_as_hdf5():
    return path.join(data_path, 'genotype_converted.hdf5')

@pytest.fixture
def genotype_imputed():
    return path.join(data_path, 'genotype_imputed.hdf5')

@pytest.fixture
def ld_folder():
    return path.join(data_path, 'ld_data')

@pytest.fixture
def genotype_as_csv():
    return path.join(data_path, 'genotype.csv')

@pytest.fixture
def pop_dataset():
    return path.join(data_path, '1k_dataset.hdf5')

@pytest.fixture
def pop_dataset_prepared():
    return path.join(data_path, '1k_dataset_unrelated.hdf5')

@pytest.fixture
def nt_map():
    return path.join(data_path, 'nt_map.pickle')
