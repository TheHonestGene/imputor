"""
Code to impute the 23andme genome for the necessary SNPs.
"""
import h5py
import scipy as sp
from itertools import *
ambig_nts = set([('A', 'T'), ('T', 'A'), ('G', 'C'), ('C', 'G')])
opp_strand_dict = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}
import gzip
import cPickle




def parse_hdf5_genotype(h5file, nt_map_file, out_h5file):
    print 'Loading NT map from file: %s'%nt_map_file
    f = open(nt_map_file, 'r')
    snp_map_dict = cPickle.load(f)
    f.close()

    print 'Parsing individual genotype: %s'%h5file
    h5f = h5py.File(h5file,'r')
    chromosomes = range(1,23) 
    
    #prepare output file
    oh5f = h5py.File(out_h5file)
    
    tot_num_parsed_snps = 0
    for chrom in chromosomes:
        print '\nWorking on chromosome %d'%chrom
        kg_chrom_str = 'chrom_%d'%chrom
        chrom_str = 'Chr%d'%chrom
        cg = h5f[chrom_str]
        sids = cg['ids'][...]
        raw_snps = cg['snps'][...]
        
        #Get the nucleotides coding map (from 1K genomes project).
        chrom_dict = snp_map_dict[kg_chrom_str]
        sid_nt_map = chrom_dict['sid_nt_map']
        n = len(sid_nt_map)
        snps = sp.repeat(-9, n)
        num_not_found = 0
        num_misunderstood = 0
        num_parsed_ok = 0
        for sid, nt in izip(sids,raw_snps):
            try:
                d = sid_nt_map[sid]
            except Exception:
                num_not_found +=1
                continue
            try:
                nt_val = d['ntm'][nt]
            except Exception:
                num_misunderstood +=1
                continue
            snps[d['snp_i']] = nt_val
            num_parsed_ok += 1
            
        print "%d SNPs weren't found and %d SNPs had unrecognizable nucleotides"%(num_not_found,num_misunderstood) 
        print "%d SNPs were parsed ok."%num_parsed_ok
        tot_num_parsed_snps +=num_parsed_ok
        #Not sure what information we need, perhaps only the SNPs?
        
        assert len(snps)==len(chrom_dict['sids'])==len(chrom_dict['positions'])==len(chrom_dict['nts']), '..bug'
        #Store information
        cg = oh5f.create_group(chrom_str)
        cg.create_dataset('snps', data=snps)
        cg.create_dataset('sids', data=chrom_dict['sids'])
        cg.create_dataset('positions', data=chrom_dict['positions'])
        cg.create_dataset('nts', data=chrom_dict['nts'])
        #Also positions..
        
        #genome_dict[chrom]={'snps':snps, } #'sids':chrom_dict['sids'], 'positions':chrom_dict['positions'], 'nts':chrom_dict['nts']}
    print 'In total %d SNPs were parsed.'%tot_num_parsed_snps
    h5f.close()
    oh5f.close()
    #return genome_dict
    

def calc_ld(ref_genotype_file, ld_matrix_h5file, window_size = 200, kgenomes_file = '/Users/bjarnivilhjalmsson/data/1Kgenomes/1K_genomes_v3_EUR_unrelated.hdf5'):
    """
    
    """
    #Load 1K genome
    kg_h5f = h5py.File(kgenomes_file)
    
    #load genotype.
    g_h5f = h5py.File(ref_genotype_file)

    mat_h5f = h5py.File(ld_matrix_h5file)

    #Figure out overlap (all genotype SNPs should be in the 1K genomes data)..
    for chrom in range(1,23):
        print 'Working on Chromosome %d'%chrom
        chrom_str1 = 'chr%d'%chrom
        kg_cg = kg_h5f[chrom_str1]
        kg_sids = kg_cg['snp_ids'][...]

        chrom_str2 = 'Chr%d'%chrom
        g_cg = g_h5f[chrom_str2]
        g_sids = g_cg['sids'][...]
        
        kg_filter = sp.in1d(kg_sids,g_sids)
        
        assert sp.sum(kg_filter)==len(g_sids), '..bug...'
        assert sp.all(kg_sids[kg_filter]==g_sids), '...bug'
        
        snps = kg_cg['snps'][...]
        snps = snps[kg_filter]
        
        snp_stds = kg_cg['snp_stds'][...]
        snp_stds = snp_stds[kg_filter]

        snp_means = kg_cg['snp_means'][...]
        snp_means = snp_means[kg_filter]
        
        norm_snps = sp.array((snps - snp_means)/snp_stds,dtype='single')
        
        #Iterate over SNPs and calculate LD
        num_snps,num_indivs = snps.shape
        ld_mat = sp.zeros((num_snps-1,window_size),dtype='single')
    
        for snp_i in range(num_snps-1):
            end_i = min(snp_i+window_size,num_snps)
            ld_mat[snp_i] = sp.dot(norm_snps[snp_i],(norm_snps[snp_i+1,end_i]).T)/num_indivs
        
        #Store things
        mat_cg = mat_h5f.create_group(chrom_str1)
        mat_cg.create_dataset('ld_mat',data=ld_mat)
    

def impute_23_and_genome(genome_dict):
    """
    """

    #For each SNP
    #  Identify a (small) window of neighboring SNPs.
    #  Identify the missing SNPs and non-missing SNPs.
    #  Invert the matrix, etc.
    #  Predict genotype
    
    
    
#For debugging purposes
if __name__=='__main__':
#     parse_hdf5_genotype('/Users/bjarnivilhjalmsson/REPOS/imputor/tests/data/test_genotype.hdf5',
#                         '/Users/bjarnivilhjalmsson/data/tmp/nt_map.pickled',
#                         '/Users/bjarnivilhjalmsson/REPOS/imputor/tests/data/test_out_genotype.hdf5')
    calc_ld('/Users/bjarnivilhjalmsson/REPOS/imputor/tests/data/test_out_genotype.hdf5', '/Users/bjarnivilhjalmsson/REPOS/imputor/tests/data/ld_mat.hdf5')