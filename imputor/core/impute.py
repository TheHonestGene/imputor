"""
Code to impute the 23andme genome for the necessary SNPs.
"""
import h5py
import scipy as sp
from itertools import *
ambig_nts = set([('A', 'T'), ('T', 'A'), ('G', 'C'), ('C', 'G')])
opp_strand_dict = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}


def coordinate_data(sum_stat_file,):
    # 1. parse nucleotide code.
    # 2. parse sum stats.
    # 3. find overlap between two.
    # 4. Store SNP key.
    pass

#Coding key
def prepare_nt_coding_key(K_genomes_snps_map,indiv_genot_file,nt_map_file):
    """
    Determines the nucleotide coding for the genotype using the 1K genome  the 10KUK
    """
    gf = h5py.File(indiv_genot_file,'r')
    kgf = h5py.File(K_genomes_snps_map,'r')
    chromosomes = range(1,23) 
    genome_dict = {}
    for chrom in chromosomes:
        kg_chrom_str = 'chrom_%d'%chrom
        chrom_str = 'Chr%d'%chrom
        
        #Get SNPs from genotype
        cg = gf[chrom_str]
        sids = cg['ids'][...]
        snps = cg['snps'][...]
        sid_dict = dict(zip(sids, snps))
        
        #Get SNP IDs from 1K genomes
        kcg = kgf[chrom_str]
        kg_sids = kcg['sids'][...]
        
        #Determine overlap between SNPs..
        kg_filter = sp.in1d(kg_sids,sids)
        kg_sids = kg_sids[kg_filter]
        kg_nts = (kcg['nts'][...])[kg_filter]
        kg_positions = (kcg['positions'][...])[kg_filter]
        
        #Check that nt are ok in genotype data, otherwise filter.
        sid_nt_map = {}
        for sid, kg_nt in izip(kg_sids, kg_nts):
            snp = sid_dict[sid]
            if snp!='--' and kg_nts not in ambig_nts:
                if snp[0]=='-' and (snp[1]== kg_nt[1] or snp[1]==kg_nt[0]):
                    sid_nt_map[sid]={'nt':kg_nt, }
        
        nts = []
        #FIXME: Nucleotide coding missing...
        
        genome_dict[chrom]={'positions':cg['positions'][...], 'ids':cg['ids'][...],'snps':snps}
    # 1. parse genotype.
    # 2. parse K genomes, create a dict with genotype, nucleotide information, etc.
    # 3. Determine and store the nt coding.



def parse_hdf5_genotype(h5file):
    h5f = h5py.File(h5file,'r')
    chromosomes = h5f.keys()
    genome_dict = {}
    for chrom in chromosomes:
        cg = h5f[chrom]
        snps = cg['snps'][...]
        nts = []
        #FIXME: Nucleotide coding missing...
        
        genome_dict[chrom]={'positions':cg['positions'][...], 'ids':cg['ids'][...],'snps':snps}
    h5f.close()
    return genome_dict
    
    
def impute_23_and_genome(genome_dict):
    """
    Liu & Stephens fast imputation of the 23andme genome.
    """
    
    #For each of the SNPs, there's a dict
    
    
#For debugging purposes
if __name__=='__main__':
    parse_hdf5_genotype()