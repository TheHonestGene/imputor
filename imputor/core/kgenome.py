"""
Methods for analysing 1000 genomes data.
"""

from itertools import *
import cPickle
import gzip
import h5py
import scipy as sp
 

__updated__ = '2016-10-12'

ambig_nts = set([('A', 'T'), ('T', 'A'), ('G', 'C'), ('C', 'G')])
opp_strand_dict = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}
Kg_nt_decoder = {1:'A', 2:'T', 3:'C', 4:'G', }

 
cloud_dir = '/Users/bjv/Dropbox/Cloud_folder/'
repos_dir = '/Users/bjv/REPOS/'

    
def gen_unrelated_eur_1k_data(input_file='/home/bjarni/TheHonestGene/faststorage/1Kgenomes/phase3/1k_genomes_hg.hdf5' ,
                              out_file='/home/bjarni/PCMA/faststorage/1_DATA/1k_genomes/1K_genomes_phase3_EUR_unrelated.hdf5',
                              maf_thres=0.01, max_relatedness=0.05, K_thinning_frac=0.1, debug=False):
    h5f = h5py.File(input_file)
    num_indivs = len(h5f['indivs']['continent'])
    eur_filter = h5f['indivs']['continent'][...] == 'EUR'
    num_eur_indivs = sp.sum(eur_filter)
    print 'Number of European individuals: %d', num_eur_indivs
    K = sp.zeros((num_eur_indivs, num_eur_indivs), dtype='single')
    num_snps = 0
    std_thres = sp.sqrt(2.0 * (1 - maf_thres) * (maf_thres))

    print 'Calculating kinship'
    for chrom in range(1, 23):
        print 'Working on Chromosome %d' % chrom
        chrom_str = 'chr%d' % chrom
        
        print 'Loading SNPs and data'
        snps = sp.array(h5f[chrom_str]['calldata']['snps'][...], dtype='int8')

        print 'Loading NTs'
        ref_nts = h5f[chrom_str]['variants']['REF'][...]
        alt_nts = h5f[chrom_str]['variants']['ALT'][...]
        
        print 'Filtering multi-allelic SNPs'
        multi_allelic_filter = sp.negative(h5f[chrom_str]['variants']['MULTI_ALLELIC'][...])
        snps = snps[multi_allelic_filter]
        ref_nts = ref_nts[multi_allelic_filter]
        alt_nts = alt_nts[multi_allelic_filter]


        if K_thinning_frac < 1:
            print 'Thinning SNPs for kinship calculation'
            thinning_filter = sp.random.random(len(snps)) < K_thinning_frac
            snps = snps[thinning_filter]
            alt_nts = alt_nts[thinning_filter]
            ref_nts = ref_nts[thinning_filter]

        print 'Filter SNPs with missing NT information'
        nt_filter = sp.in1d(ref_nts, ok_nts)
        nt_filter = nt_filter * sp.in1d(alt_nts, ok_nts)
        if sp.sum(nt_filter) < len(nt_filter):
            snps = snps[nt_filter]

        print 'Filtering non-European individuals'
        snps = snps[:, eur_filter]

        print 'Filtering SNPs with MAF <', maf_thres
        snp_stds = sp.std(snps, 1)
        maf_filter = snp_stds.flatten() > std_thres
        snps = snps[maf_filter]
        snp_stds = snp_stds[maf_filter]
        
        print '%d SNPs remaining after all filtering steps.' % len(snps)

        print 'Normalizing SNPs'
        snp_means = sp.mean(snps, 1)
        norm_snps = (snps - snp_means[sp.newaxis].T) / snp_stds[sp.newaxis].T
        
        print 'Updating kinship'        
        K += sp.dot(norm_snps.T, norm_snps)
        num_snps += len(norm_snps)
        assert sp.isclose(sp.sum(sp.diag(K)) / (num_snps * num_eur_indivs), 1.0)

    K = K / float(num_snps)
    print 'Kinship calculation done using %d SNPs\n' % num_snps
    
    # Filter individuals
    print 'Filtering individuals'
    keep_indiv_set = set(range(num_eur_indivs))
    for i in range(num_eur_indivs):
        if i in keep_indiv_set:
            for j in range(i + 1, num_eur_indivs):
                if K[i, j] > max_relatedness:
                    if j in keep_indiv_set:
                        keep_indiv_set.remove(j)
    keep_indivs = list(keep_indiv_set)
    keep_indivs.sort()
    print 'Retained %d individuals\n' % len(keep_indivs)
    
    # Checking that everything is ok!
    K_ok = K[keep_indivs]
    K_ok = K_ok[:, keep_indivs]
    assert (K_ok - sp.tril(K_ok)).max() < max_relatedness

    indiv_filter = sp.zeros(num_indivs, dtype='bool8')
    indiv_filter[(sp.arange(num_indivs)[eur_filter])[keep_indivs]] = 1
    
    assert sp.sum(indiv_filter) == len(keep_indivs)
    
    # Store in new file
    print 'Now storing data.'
    oh5f = h5py.File(out_file, 'w')
    indiv_ids = h5f['indivs']['indiv_ids'][indiv_filter]
    oh5f.create_dataset('indiv_ids', data=indiv_ids)    
    for chrom in range(1, 23):
        print 'Working on Chromosome %d' % chrom
        chrom_str = 'chr%d' % chrom
        
        print 'Loading SNPs and data'
        snps = sp.array(h5f[chrom_str]['calldata']['snps'][...], dtype='int8')
        snp_ids = h5f[chrom_str]['variants']['ID'][...]
        positions = h5f[chrom_str]['variants']['POS'][...]

        print 'Loading NTs'
        ref_nts = h5f[chrom_str]['variants']['REF'][...]
        alt_nts = h5f[chrom_str]['variants']['ALT'][...]
        
        print 'Filtering multi-allelic SNPs'
        multi_allelic_filter = sp.negative(h5f[chrom_str]['variants']['MULTI_ALLELIC'][...])
        snps = snps[multi_allelic_filter]
        ref_nts = ref_nts[multi_allelic_filter]
        alt_nts = alt_nts[multi_allelic_filter]
        positions = positions[multi_allelic_filter]
        snp_ids = snp_ids[multi_allelic_filter]

        print 'Filter individuals'
        snps = snps[:, indiv_filter]
        
        print 'Filter SNPs with missing NT information'
        nt_filter = sp.in1d(ref_nts, ok_nts)
        nt_filter = nt_filter * sp.in1d(alt_nts, ok_nts)
        if sp.sum(nt_filter) < len(nt_filter):
            snps = snps[nt_filter]
            ref_nts = ref_nts[nt_filter]
            alt_nts = alt_nts[nt_filter]
            positions = positions[nt_filter]
            snp_ids = snp_ids[nt_filter]
        
        print 'filter monomorphic SNPs'
        snp_stds = sp.std(snps, 1)
        mono_morph_filter = snp_stds > 0
        snps = snps[mono_morph_filter]
        ref_nts = ref_nts[mono_morph_filter]
        alt_nts = alt_nts[mono_morph_filter]
        positions = positions[mono_morph_filter]
        snp_ids = snp_ids[mono_morph_filter]
        snp_stds = snp_stds[mono_morph_filter]

        snp_means = sp.mean(snps, 1)

        if debug:
            if K_thinning_frac < 1:
                print 'Thinning SNPs for kinship calculation'
                thinning_filter = sp.random.random(len(snps)) < K_thinning_frac
                k_snps = snps[thinning_filter]
                k_snp_stds = snp_stds[thinning_filter]

    
            print 'Filtering SNPs with MAF <', maf_thres
            maf_filter = k_snp_stds.flatten() > std_thres
            k_snps = k_snps[maf_filter]
            k_snp_stds = k_snp_stds[maf_filter]
            k_snp_means = sp.mean(k_snps)

            print 'Verifying that the Kinship makes sense'
            norm_snps = (k_snps - k_snp_means[sp.newaxis].T) / k_snp_stds[sp.newaxis].T
            K = sp.dot(norm_snps.T, norm_snps)
            num_snps += len(norm_snps)
            if sp.isclose(sp.sum(sp.diag(K)) / (num_snps * num_eur_indivs), 1.0) and (K - sp.tril(K)).max() < (max_relatedness * 1.5):
                print 'It looks OK!'
            else:
                raise Exception('Kinship looks wrong?')
        

        nts = sp.array([[nt1, nt2] for nt1, nt2 in izip(ref_nts, alt_nts)])

        print 'Writing to disk'
        cg = oh5f.create_group(chrom_str)
        cg.create_dataset('snps', data=snps)
        cg.create_dataset('snp_means', data=snp_means[sp.newaxis].T)
        cg.create_dataset('snp_stds', data=snp_stds[sp.newaxis].T)
        cg.create_dataset('snp_ids', data=snp_ids)
        cg.create_dataset('positions', data=positions)
        cg.create_dataset('nts', data=nts)
        oh5f.flush()
        print 'Done writing to disk'
        
#         centimorgans = h5f[chrom_str]['centimorgans'][...]
#         cg.create_dataset('centimorgans',data=centimorgans)
#         
#         centimorgan_rates = h5f[chrom_str]['centimorgan_rates'][...]
#         cg.create_dataset('centimorgan_rates',data=centimorgan_rates)
        
    oh5f.close()
    h5f.close()
    print 'Done'
    
   
def gen_1k_test_genotypes(kg_file=cloud_dir + 'Data/1Kgenomes/1K_genomes_v3_EUR_unrelated2.hdf5',
                          nt_map_file=cloud_dir + 'tmp/nt_map.pickled', out_prefix=cloud_dir + 'tmp/1k_ind'):
    """
    Generates 1K genotypes in the internal genotype format for validation and other purposes.
    """
    print 'Loading NT map from file: %s' % nt_map_file
    f = open(nt_map_file, 'r')
    snp_map_dict = cPickle.load(f)
    f.close()
    
    h5f = h5py.File(kg_file)
    kg_indivs = h5f['indiv_ids'][...]
    chromosomes = range(1, 23) 
    
    
    # Figure out 1K SNP filter
    chrom_filter_dict = {}
    
    for chrom in chromosomes:
        kg_chrom_str = 'chr%d' % chrom
        cg = h5f[kg_chrom_str]
        sids = cg['snp_ids'][...]
        
        # Get the nucleotides coding map (from 1K genomes project).
        chrom_dict = snp_map_dict[kg_chrom_str]
        ok_sids = chrom_dict['sids']
        snps_filter = sp.in1d(sids, ok_sids)
        chrom_filter_dict[chrom] = snps_filter
        
        reorder_snps = False
        filtered_sids = sids[snps_filter]
        assert len(filtered_sids) == len(ok_sids), '.... bug'
        if not sp.all(filtered_sids == ok_sids):
            sids_indices_dict = {}
            for i, sid in enumerate(sids):
                sids_indices_dict[sid] = i
            snp_order = []
            for sid in ok_sids:
                snp_order.append(sids_indices_dict[sid])
            filtered_sids = filtered_sids[snp_order]
            assert sp.all(filtered_sids == ok_sids), '... bug'
            reorder_snps = True
    
    
    for ind_i in range(10):  # len(kg_indivs)):
        print 'Generating genotype for individual: %d ' % ind_i
        
        # prepare output file
        out_h5fn = out_prefix + '_%d.hdf5' % ind_i
        out_h5f = h5py.File(out_h5fn)
        
        for chrom in chromosomes:
#             print '\nWorking on chromosome %d'%chrom
            kg_chrom_str = 'chr%d' % chrom
            chrom_str = 'Chr%d' % chrom
            cg = h5f[kg_chrom_str]
            
            snps_filter = chrom_filter_dict[chrom]

            sids = (cg['snp_ids'][...])[snps_filter]
            snps = (cg['snps'][...])[:, ind_i]
            snps = snps[snps_filter]
            positions = (cg['positions'][...])[snps_filter]
            nts = (cg['nts'][...])[snps_filter]
            
            if reorder_snps:
                snps = snps[snp_order]
                positions = positions[snp_order]
                nts = nts[snp_order]
                sids = sids[snp_order]
            
            assert len(snps) == len(sids) == len(positions) == len(nts), '..bug'
            # Store information
            ocg = out_h5f.create_group  (chrom_str)
            ocg.create_dataset('snps', data=snps)
            ocg.create_dataset('sids', data=sids)
            ocg.create_dataset('positions', data=positions)
            ocg.create_dataset('nts', data=nts)

        out_h5f.close()
    h5f.close()
        # return genome_dict        



# Coding key
# def prepare_nt_coding_key(K_genomes_snps_map, indiv_genot_file, nt_map_file):
#     """
#     Determines the nucleotide coding for the genotype using the 1K genome  the 10KUK
#     """
#     gf = h5py.File(indiv_genot_file,'r')
#     kgf = h5py.File(K_genomes_snps_map,'r')
#     chromosomes = range(1,23) 
#     snp_map_dict = {}
#     num_snps = 0
#     for chrom in chromosomes:
#         print 'Working on chromosome %d'%chrom
#         kg_chrom_str = 'chr%d'%chrom
#         chrom_str = 'Chr%d'%chrom
#         
#         #Get SNPs from genotype
#         cg = gf[chrom_str]
#         sids = cg['ids'][...]
#         snps = cg['snps'][...]
#         sid_dict = dict(zip(sids, snps))
#         
#         #Get SNP IDs from 1K genomes
#         kcg = kgf[kg_chrom_str]
#         kg_sids = kcg['snp_ids'][...]
#         
#         #Determine overlap between SNPs..
#         kg_filter = sp.in1d(kg_sids,sids)
#         kg_sids = kg_sids[kg_filter]
#         kg_nts = (kcg['nts'][...])[kg_filter]
#         kg_positions = (kcg['positions'][...])[kg_filter]
#         
#         #Check that nt are ok in genotype data, otherwise filter.
#         sid_nt_map = {}
#         positions = []
#         ok_sids = []
#         nts = []
#         snp_i = 0
#         for sid, kg_nt, kg_pos in izip(kg_sids, kg_nts, kg_positions):
#             snp = sid_dict[sid]
#             if tuple(kg_nt) not in ambig_nts:
#                 # All possible (and allowed) nucleotides strings 
#                 ntm = {}
#                 ntm['--']=-9
#                 ntm['-'+kg_nt[0]]=-9
#                 ntm['-'+kg_nt[1]]=-9
#                 ntm[kg_nt[0]+'-']=-9
#                 ntm[kg_nt[1]+'-']=-9
#                 ntm[kg_nt[0]+kg_nt[0]]=0
#                 ntm[kg_nt[1]+kg_nt[0]]=1
#                 ntm[kg_nt[0]+kg_nt[1]]=1
#                 ntm[kg_nt[1]+kg_nt[1]]=2
#                 sid_nt_map[sid]={'ntm':ntm, 'snp_i':snp_i}
#                 positions.append(kg_pos)
#                 nts.append(kg_nt)
#                 ok_sids.append(sid)
#                 snp_i += 1
#         
#         num_snps += len(sid_nt_map)
#         snp_map_dict[kg_chrom_str]={'sid_nt_map':sid_nt_map, 'positions':positions, 'nts':nts, 'sids':ok_sids}
#     
#     print 'Found %d SNPs'%num_snps
#     print 'Writing to file'
#     f = open(nt_map_file, 'wb')
#     cPickle.dump(snp_map_dict, f, protocol=2)
#     f.close()
#     return snp_map_dict
        
# For debugging purposes
if __name__ == '__main__':        
    gen_1k_test_genotypes()


