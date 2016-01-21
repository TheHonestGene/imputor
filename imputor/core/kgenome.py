import h5py
import scipy as sp
from itertools import *
ambig_nts = set([('A', 'T'), ('T', 'A'), ('G', 'C'), ('C', 'G')])
opp_strand_dict = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}
import gzip
import cPickle

    
def gen_unrelated_eur_1k_data(out_file='/Users/bjarnivilhjalmsson/data/1Kgenomes/1K_genomes_v3_EUR_unrelated.hdf5'):
    h5f = h5py.File('/Users/bjarnivilhjalmsson/data/1Kgenomes/1K_genomes_v3.hdf5')
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

    #Map prefix
    KGenomes_prefix='/Users/bjarnivilhjalmsson/data/1Kgenomes/', 
    
    #Filter ambiguous SNPs and store in new file
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
        snps = snps[:,keep_indivs]
        snp_stds = sp.std(snps,1)
        mono_morph_filter = snp_stds>0
        snps = snps[mono_morph_filter]
        
        #Load map and double check..
        print 'Load  SNP map to get NTs and filter ambiguous SNPs'
        fn = '%sALL_1000G_phase1integrated_v3_chr%d_impute.legend.gz'%(KGenomes_prefix,chrom_i)
        ok_sids = set()
        sids_nt_map = {}
        with gzip.open(fn) as f:
            f.next()
            line_i=0
            for line in f:
                l = line.split()
                nt1=l[2]
                nt2=l[3]
                if nt1 not in valid_nts:
                    continue
                if nt2 not in valid_nts:
                    continue
                if filter_ambiguous and (nt1,nt2) in ambig_nts:
                    continue
                line_i +=1
                ok_sids.add(l[0])
                sids_nt_map[l[0]]=(nt1,nt2)
        
        #FINISH!!!!

        
        cg = oh5f.create_group(chrom_str)
        
        cg.create_dataset('snps',data=snps)

        snp_stds = snp_stds[mono_morph_filter]
        snp_means = sp.mean(snps,1)
        cg.create_dataset('snp_means',data=snp_means[sp.newaxis].T)
        cg.create_dataset('snp_stds',data=snp_stds[sp.newaxis].T)
       
        snp_ids = h5f[chrom_str]['snp_ids'][...]
        snp_ids = snp_ids[mono_morph_filter]
        cg.create_dataset('snp_ids',data=snp_ids)
        
        positions = h5f[chrom_str]['positions'][...]
        positions = positions[mono_morph_filter]
        cg.create_dataset('positions',data=positions)
        
        eur_maf = h5f[chrom_str]['eur_maf'][...]
        eur_maf = eur_maf[mono_morph_filter]
        cg.create_dataset('eur_maf',data=eur_maf)
        
        nts = h5f[chrom_str]['nts'][...]
        nts = nts[mono_morph_filter]
        cg.create_dataset('nts',data=nts)
        
        centimorgans = h5f[chrom_str]['centimorgans'][...]
        cg.create_dataset('centimorgans',data=centimorgans)
        
        centimorgan_rates = h5f[chrom_str]['centimorgan_rates'][...]
        cg.create_dataset('centimorgan_rates',data=centimorgan_rates)
        
    oh5f.close()
    h5f.close()
   




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
        
#For debugging purposes
if __name__=='__main__':        
    prepare_nt_coding_key('/Users/bjarnivilhjalmsson/data/1Kgenomes/1K_genomes_v3_EUR_unrelated.hdf5',
                          '/Users/bjarnivilhjalmsson/REPOS/imputor/tests/data/test_genotype.hdf5',
                          '/Users/bjarnivilhjalmsson/data/tmp/nt_map2.pickled')


