"""
Code to impute the 23andme genome for the necessary SNPs.
"""
import h5py
import scipy as sp
from itertools import *
ambig_nts = set([('A', 'T'), ('T', 'A'), ('G', 'C'), ('C', 'G')])
opp_strand_dict = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}
Kg_nt_decoder = {1:'A', 2:'T', 3:'C', 4:'G', }

import gzip
import cPickle


#Coding key
def prepare_nt_coding_key(K_genomes_snps_map, indiv_genot_file, nt_map_file):
    """
    Determines the nucleotide coding for the genotype using the 1K genome  
    """
    gf = h5py.File(indiv_genot_file,'r')
    kgf = h5py.File(K_genomes_snps_map,'r')
    chromosomes = range(1,23) 
    snp_map_dict = {}
    num_snps = 0
    for chrom in chromosomes:
        print 'Working on chromosome %d'%chrom
#         kg_chrom_str = 'chrom_%d'%chrom
        kg_chrom_str = 'chr%d'%chrom
        chrom_str = 'Chr%d'%chrom
        
        #Get SNPs from genotype
        cg = gf[chrom_str]
        sids = cg['ids'][...]
        snps = cg['snps'][...]
        sid_dict = dict(zip(sids, snps))
        
        #Get SNP IDs from 1K genomes
        kcg = kgf[kg_chrom_str]
        kg_sids = kcg['snp_ids'][...]
        
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
        kg_chrom_str = 'chr%d'%chrom
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
    
    
def gen_unrelated_eur_1k_data(out_file='/Users/bjv/Dropbox/Cloud_folder/Data/1Kgenomes/1K_genomes_v3_EUR_unrelated2.hdf5'):
    h5f = h5py.File('/Users/bjv/Dropbox/Cloud_folder/Data/1Kgenomes/1K_genomes_v3.hdf5')
    eur_filter = h5f['indivs']['continent'][...]=='EUR'
    num_indivs = sp.sum(eur_filter)
    K = sp.zeros((num_indivs,num_indivs), dtype='single')
    num_snps = 0
    print 'Calculating kinship'
    for chrom in range(1,23):
        print 'Working on Chromosome %d'%chrom
        chrom_str = 'chr%d'%chrom
        print 'Loading SNPs'
        snps = h5f[chrom_str]['raw_snps'][...]
        #filter non-europeans.
        print 'Filtering non-European individuals'
        snps = snps[:,eur_filter]
        print 'Filtering monomorphic SNPs'
        snp_stds = sp.std(snps,1)
        mono_morph_filter = snp_stds>0
        snps = snps[mono_morph_filter]
        snp_stds = snp_stds[mono_morph_filter]
        print 'Filter SNPs with missing NT information'
        nts = h5f[chrom_str]['nts'][...]
        nts = nts[mono_morph_filter]
        nt_filter = sp.all(nts>0,1)
        snps = snps[nt_filter]
        snp_stds = snp_stds[nt_filter]
        print 'Normalizing SNPs'
        snp_means = sp.mean(snps,1)
        norm_snps = (snps - snp_means[sp.newaxis].T)/snp_stds[sp.newaxis].T
        print 'Updating kinship'        
        K += sp.dot(norm_snps.T,norm_snps)
        num_snps += len(norm_snps)
    K = K/float(num_snps)
    print 'Kinship calculation done using %d SNPs\n'%num_snps
    
    
    #Filter individuals
    print 'Filtering individuals'
    keep_indiv_set = set(range(num_indivs))
    for i in range(num_indivs):
        if i in keep_indiv_set:
            for j in range(i+1,num_indivs):
                if K[i,j]>0.05:
                    if j in keep_indiv_set:
                        keep_indiv_set.remove(j)
    keep_indivs = list(keep_indiv_set)
    keep_indivs.sort()
    print 'Retained %d individuals\n'%len(keep_indivs)

    #Store in new file
    print 'Now storing data.'
    oh5f = h5py.File(out_file,'w')
    indiv_ids = h5f['indivs']['indiv_ids'][eur_filter]
    indiv_ids = indiv_ids[keep_indivs]
    oh5f.create_dataset('indiv_ids',data =indiv_ids)    
    for chrom in range(1,23):
        print 'Working on Chromosome %d'%chrom
        chrom_str = 'chr%d'%chrom
        snps = h5f[chrom_str]['raw_snps'][...]
        #filter non-europeans.
        snps = snps[:,eur_filter]
        #Filter related
        snps = snps[:,keep_indivs]
        #filter monomorphic SNPs
        snp_stds = sp.std(snps,1)
        mono_morph_filter = snp_stds>0
        snps = snps[mono_morph_filter]
        #filter SNPs w missing NT values
        nts = h5f[chrom_str]['nts'][...]
        nts = nts[mono_morph_filter]
        nt_filter = sp.all(nts>0,1)
        nts = nts[nt_filter]
        snps = snps[nt_filter]
        
        cg = oh5f.create_group(chrom_str)
        cg.create_dataset('snps',data=snps)

        snp_stds = snp_stds[mono_morph_filter]
        snp_stds = snp_stds[nt_filter]
        snp_means = sp.mean(snps,1)
        cg.create_dataset('snp_means',data=snp_means[sp.newaxis].T)
        cg.create_dataset('snp_stds',data=snp_stds[sp.newaxis].T)
       
        snp_ids = h5f[chrom_str]['snp_ids'][...]
        snp_ids = snp_ids[mono_morph_filter]
        snp_ids = snp_ids[nt_filter]
        cg.create_dataset('snp_ids',data=snp_ids)
        
        positions = h5f[chrom_str]['positions'][...]
        positions = positions[mono_morph_filter]
        positions = positions[nt_filter]
        cg.create_dataset('positions',data=positions)
        
        eur_maf = h5f[chrom_str]['eur_maf'][...]
        eur_maf = eur_maf[mono_morph_filter]
        eur_maf = eur_maf[nt_filter]
        cg.create_dataset('eur_maf',data=eur_maf)
        
        decoded_nts = []
        for nt in nts:
            nt1 = nt[0]
            nt2 = nt[1]
            decoded_nts.append([Kg_nt_decoder[nt1],Kg_nt_decoder[nt2]])
        cg.create_dataset('nts',data=decoded_nts)
        
        centimorgans = h5f[chrom_str]['centimorgans'][...]
        cg.create_dataset('centimorgans',data=centimorgans)
        
        centimorgan_rates = h5f[chrom_str]['centimorgan_rates'][...]
        cg.create_dataset('centimorgan_rates',data=centimorgan_rates)
        
    oh5f.close()
    h5f.close()
   
    
def calc_ld(ref_genotype_file, ld_matrix_h5file, window_size = 200, kgenomes_file = '/Users/bjv/Dropbox/Cloud_folder/Data/1Kgenomes/1K_genomes_v3_EUR_unrelated2.hdf5'):
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
            #Deal with boundaries..
            if snp_i+window_size>num_snps:
                
            ld_mat[snp_i] = sp.dot(norm_snps[snp_i],(norm_snps[snp_i+1:end_i+1]).T)/num_indivs
        
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
    pass
    

    
#For debugging purposes
if __name__=='__main__':
#     Filter related indivs
#     gen_unrelated_eur_1k_data()
    
#     prepare_nt_coding_key('/Users/bjv/Dropbox/Cloud_folder/Data/1Kgenomes/1K_genomes_v3_EUR_unrelated2.hdf5',
#                           '/Users/bjv/REPOS/imputor/tests/data/test_genotype.hdf5',
#                           '/Users/bjv/Dropbox/Cloud_folder/tmp/nt_map.pickled')
#     parse_hdf5_genotype('/Users/bjv/REPOS/imputor/tests/data/test_genotype.hdf5',
#                          '/Users/bjv/Dropbox/Cloud_folder/tmp/nt_map.pickled',
#                          '/Users/bjv/REPOS/imputor/tests/data/test_out_genotype.hdf5')
    calc_ld('/Users/bjv/REPOS/imputor/tests/data/test_out_genotype.hdf5', '/Users/bjv/REPOS/imputor/tests/data/ld_mat.hdf5')
    
    
