"""
Methods for analysing 1000 genomes data.
"""

import h5py
import scipy as sp
import logging

try:
    from itertools import izip as zip
except ImportError:  # will be 3.x series
    pass

try:
    import cPickle
except ImportError:
    import pickle as cPickle

ambig_nts = set([('A', 'T'), ('T', 'A'), ('G', 'C'), ('C', 'G')])
opp_strand_dict = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}
Kg_nt_decoder = {1:'A', 2:'T', 3:'C', 4:'G', }

ok_nts = ['A', 'C', 'G', 'T']

log = logging.getLogger(__name__)

cloud_dir = '/Users/bjv/Dropbox/Cloud_folder/'
repos_dir = '/Users/bjv/REPOS/'


def gen_unrelated_eur_1k_data(input_file='/home/bjarni/TheHonestGene/faststorage/1Kgenomes/phase3/1k_genomes_hg.hdf5' ,
                              out_file='/home/bjarni/PCMA/faststorage/1_DATA/1k_genomes/1K_genomes_phase3_EUR_unrelated.hdf5',
                              maf_thres=0.01, max_relatedness=0.05, K_thinning_frac=0.1, debug=False):
    h5f = h5py.File(input_file)
    num_indivs = len(h5f['indivs']['continent'])
    eur_filter = h5f['indivs']['continent'][...] == 'EUR'
    num_eur_indivs = sp.sum(eur_filter)
    log.info('Number of European individuals: %d', num_eur_indivs)
    K = sp.zeros((num_eur_indivs, num_eur_indivs), dtype='single')
    num_snps = 0
    std_thres = sp.sqrt(2.0 * (1 - maf_thres) * (maf_thres))

    log.info('Calculating kinship')
    for chrom in range(1, 23):
        log.info('Working on Chromosome %d' % chrom)
        chrom_str = 'chr%d' % chrom

        log.info('Loading SNPs and data')
        snps = sp.array(h5f[chrom_str]['calldata']['snps'][...], dtype='int8')

        log.info('Loading NTs')
        ref_nts = h5f[chrom_str]['variants']['REF'][...]
        alt_nts = h5f[chrom_str]['variants']['ALT'][...]

        log.info('Filtering multi-allelic SNPs')
        multi_allelic_filter = sp.negative(h5f[chrom_str]['variants']['MULTI_ALLELIC'][...])
        snps = snps[multi_allelic_filter]
        ref_nts = ref_nts[multi_allelic_filter]
        alt_nts = alt_nts[multi_allelic_filter]


        if K_thinning_frac < 1:
            log.info('Thinning SNPs for kinship calculation')
            thinning_filter = sp.random.random(len(snps)) < K_thinning_frac
            snps = snps[thinning_filter]
            alt_nts = alt_nts[thinning_filter]
            ref_nts = ref_nts[thinning_filter]

        log.info('Filter SNPs with missing NT information')
        nt_filter = sp.in1d(ref_nts, ok_nts)
        nt_filter = nt_filter * sp.in1d(alt_nts, ok_nts)
        if sp.sum(nt_filter) < len(nt_filter):
            snps = snps[nt_filter]

        log.info('Filtering non-European individuals')
        snps = snps[:, eur_filter]

        log.info('Filtering SNPs with MAF <', maf_thres)
        snp_stds = sp.std(snps, 1)
        maf_filter = snp_stds.flatten() > std_thres
        snps = snps[maf_filter]
        snp_stds = snp_stds[maf_filter]

        log.info('%d SNPs remaining after all filtering steps.' % len(snps))

        log.info('Normalizing SNPs')
        snp_means = sp.mean(snps, 1)
        norm_snps = (snps - snp_means[sp.newaxis].T) / snp_stds[sp.newaxis].T

        log.info('Updating kinship')
        K += sp.dot(norm_snps.T, norm_snps)
        num_snps += len(norm_snps)
        assert sp.isclose(sp.sum(sp.diag(K)) / (num_snps * num_eur_indivs), 1.0)

    K = K / float(num_snps)
    log.info('Kinship calculation done using %d SNPs\n' % num_snps)

    # Filter individuals
    log.info('Filtering individuals')
    keep_indiv_set = set(range(num_eur_indivs))
    for i in range(num_eur_indivs):
        if i in keep_indiv_set:
            for j in range(i + 1, num_eur_indivs):
                if K[i, j] > max_relatedness:
                    if j in keep_indiv_set:
                        keep_indiv_set.remove(j)
    keep_indivs = list(keep_indiv_set)
    keep_indivs.sort()
    log.info('Retained %d individuals\n' % len(keep_indivs))

    # Checking that everything is ok!
    K_ok = K[keep_indivs]
    K_ok = K_ok[:, keep_indivs]
    assert (K_ok - sp.tril(K_ok)).max() < max_relatedness

    indiv_filter = sp.zeros(num_indivs, dtype='bool8')
    indiv_filter[(sp.arange(num_indivs)[eur_filter])[keep_indivs]] = 1

    assert sp.sum(indiv_filter) == len(keep_indivs)

    # Store in new file
    log.info('Now storing data.')
    oh5f = h5py.File(out_file, 'w')
    indiv_ids = h5f['indivs']['indiv_ids'][indiv_filter]
    oh5f.create_dataset('indiv_ids', data=indiv_ids)
    for chrom in range(1, 23):
        log.info('Working on Chromosome %d' % chrom)
        chrom_str = 'chr%d' % chrom

        log.info('Loading SNPs and data')
        snps = sp.array(h5f[chrom_str]['calldata']['snps'][...], dtype='int8')
        snp_ids = h5f[chrom_str]['variants']['ID'][...]
        positions = h5f[chrom_str]['variants']['POS'][...]

        log.info('Loading NTs')
        ref_nts = h5f[chrom_str]['variants']['REF'][...]
        alt_nts = h5f[chrom_str]['variants']['ALT'][...]

        log.info('Filtering multi-allelic SNPs')
        multi_allelic_filter = sp.negative(h5f[chrom_str]['variants']['MULTI_ALLELIC'][...])
        snps = snps[multi_allelic_filter]
        ref_nts = ref_nts[multi_allelic_filter]
        alt_nts = alt_nts[multi_allelic_filter]
        positions = positions[multi_allelic_filter]
        snp_ids = snp_ids[multi_allelic_filter]

        log.info('Filter individuals')
        snps = snps[:, indiv_filter]

        log.info('Filter SNPs with missing NT information')
        nt_filter = sp.in1d(ref_nts, ok_nts)
        nt_filter = nt_filter * sp.in1d(alt_nts, ok_nts)
        if sp.sum(nt_filter) < len(nt_filter):
            snps = snps[nt_filter]
            ref_nts = ref_nts[nt_filter]
            alt_nts = alt_nts[nt_filter]
            positions = positions[nt_filter]
            snp_ids = snp_ids[nt_filter]

        log.info('filter monomorphic SNPs')
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
                log.info('Thinning SNPs for kinship calculation')
                thinning_filter = sp.random.random(len(snps)) < K_thinning_frac
                k_snps = snps[thinning_filter]
                k_snp_stds = snp_stds[thinning_filter]


            log.info('Filtering SNPs with MAF <', maf_thres)
            maf_filter = k_snp_stds.flatten() > std_thres
            k_snps = k_snps[maf_filter]
            k_snp_stds = k_snp_stds[maf_filter]
            k_snp_means = sp.mean(k_snps)

            log.info('Verifying that the Kinship makes sense')
            norm_snps = (k_snps - k_snp_means[sp.newaxis].T) / k_snp_stds[sp.newaxis].T
            K = sp.dot(norm_snps.T, norm_snps)
            num_snps += len(norm_snps)
            if sp.isclose(sp.sum(sp.diag(K)) / (num_snps * num_eur_indivs), 1.0) and (K - sp.tril(K)).max() < (max_relatedness * 1.5):
                log.info('It looks OK!')
            else:
                raise Exception('Kinship looks wrong?')


        nts = sp.array([[nt1, nt2] for nt1, nt2 in zip(ref_nts, alt_nts)])

        log.info('Writing to disk')
        cg = oh5f.create_group(chrom_str)
        cg.create_dataset('snps', data=snps)
        cg.create_dataset('snp_means', data=snp_means[sp.newaxis].T)
        cg.create_dataset('snp_stds', data=snp_stds[sp.newaxis].T)
        cg.create_dataset('snp_ids', data=snp_ids)
        cg.create_dataset('positions', data=positions)
        cg.create_dataset('nts', data=nts)
        oh5f.flush()
        log.info('Done writing to disk')


    oh5f.close()
    h5f.close()
    log.info('Done')


def gen_1k_test_genotypes(kg_file=cloud_dir + 'Data/1Kgenomes/1K_genomes_v3_EUR_unrelated2.hdf5',
                          nt_map_file=cloud_dir + 'tmp/nt_map.pickled', out_prefix=cloud_dir + 'tmp/1k_ind'):
    """
    Generates 1K genotypes in the internal genotype format for validation and other purposes.
    """
    log.info('Loading NT map from file: %s' % nt_map_file)
    f = open(nt_map_file, 'r')
    snp_map_dict = cPickle.load(f)
    f.close()

    h5f = h5py.File(kg_file)
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
        log.info('Generating genotype for individual: %d ' % ind_i)

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


# For debugging purposes
if __name__ == '__main__':
    gen_1k_test_genotypes()


