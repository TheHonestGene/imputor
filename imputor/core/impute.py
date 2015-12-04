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



#Coding key
def prepare_nt_coding_key(K_genomes_snps_map, indiv_genot_file, nt_map_file):
    """
    Determines the nucleotide coding for the genotype using the 1K genome  the 10KUK
    """
    gf = h5py.File(indiv_genot_file,'r')
    kgf = h5py.File(K_genomes_snps_map,'r')
    chromosomes = range(1,23) 
    snp_map_dict = {}
    num_snps = 0
    for chrom in chromosomes:
        print 'Working on chromosome %d'%chrom
        kg_chrom_str = 'chrom_%d'%chrom
        chrom_str = 'Chr%d'%chrom
        
        #Get SNPs from genotype
        cg = gf[chrom_str]
        sids = cg['ids'][...]
        snps = cg['snps'][...]
        sid_dict = dict(zip(sids, snps))
        
        #Get SNP IDs from 1K genomes
        kcg = kgf[kg_chrom_str]
        kg_sids = kcg['sids'][...]
        
        #Determine overlap between SNPs..
        kg_filter = sp.in1d(kg_sids,sids)
        kg_sids = kg_sids[kg_filter]
        kg_nts = (kcg['nts'][...])[kg_filter]
        kg_positions = (kcg['positions'][...])[kg_filter]
        
        #Check that nt are ok in genotype data, otherwise filter.
        sid_nt_map = {}
        positions = []
        ok_sids = []
        nts = []
        snp_i = 0
        for sid, kg_nt, kg_pos in izip(kg_sids, kg_nts, kg_positions):
            snp = sid_dict[sid]
            if tuple(kg_nt) not in ambig_nts:
                # All possible (and allowed) nucleotides strings 
                ntm = {}
                ntm['--']=-9
                ntm['-'+kg_nt[0]]=-9
                ntm['-'+kg_nt[1]]=-9
                ntm[kg_nt[0]+'-']=-9
                ntm[kg_nt[1]+'-']=-9
                ntm[kg_nt[0]+kg_nt[0]]=0
                ntm[kg_nt[1]+kg_nt[0]]=1
                ntm[kg_nt[0]+kg_nt[1]]=1
                ntm[kg_nt[1]+kg_nt[1]]=2
                sid_nt_map[sid]={'ntm':ntm, 'snp_i':snp_i}
                positions.append(kg_pos)
                nts.append(kg_nt)
                ok_sids.append(sid)
                snp_i += 1
        
        num_snps += len(sid_nt_map)
        snp_map_dict[kg_chrom_str]={'sid_nt_map':sid_nt_map, 'positions':positions, 'nts':nts, 'sids':ok_sids}
    
    print 'Found %d SNPs'%num_snps
    print 'Writing to file'
    f = open(nt_map_file, 'wb')
    cPickle.dump(snp_map_dict, f, protocol=2)
    f.close()
    return snp_map_dict
        

def parse_hdf5_genotype(h5file, nt_map_file):
    print 'Loading NT map from file: %s'%nt_map_file
    f = open(nt_map_file, 'r')
    snp_map_dict = cPickle.load(f)
    f.close()

    print 'Parsing individual genotype: %s'%h5file
    h5f = h5py.File(h5file,'r')
    chromosomes = range(1,23) 
    genome_dict = {}
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
        genome_dict[chrom]={'snps':snps, } #'sids':chrom_dict['sids'], 'positions':chrom_dict['positions'], 'nts':chrom_dict['nts']}
    print 'In total %d SNPs were parsed.'%tot_num_parsed_snps
    h5f.close()
    return genome_dict
    
    
def impute_23_and_genome(genome_dict):
    """
    Liu & Stephens fast imputation of the 23andme genome.
    """
    
    #For each of the SNPs, there's a dict
    
    
#For debugging purposes
if __name__=='__main__':
#     prepare_nt_coding_key('/Users/bjarnivilhjalmsson/data/1Kgenomes/snps.hdf5',
#                           '/Users/bjarnivilhjalmsson/REPOS/imputor/tests/data/test_genotype.hdf5',
#                           '/Users/bjarnivilhjalmsson/data/tmp/nt_map.pickled')
    parse_hdf5_genotype('/Users/bjarnivilhjalmsson/REPOS/imputor/tests/data/test_genotype.hdf5',
                        '/Users/bjarnivilhjalmsson/data/tmp/nt_map.pickled')