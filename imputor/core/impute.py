"""
Code to impute the 23andme genome for the necessary SNPs.

"""
from functools import reduce  # # py3 does not have it 
from sys import version_info
import gzip
import logging
import operator
import random

from scipy import linalg 
import h5py

import scipy as sp


try:
    from itertools import izip as zip
except ImportError:  # will be 3.x series
    pass
try:
    import cPickle as pickle
except ImportError:  # will be 3.series
    import pickle



ambig_nts = set([(b'A', b'T'), (b'T', b'A'), (b'G', b'C'), (b'C', b'G')])

log = logging.getLogger(__name__)

def _get_chunk_length(data):
    rowsize = data.dtype.itemsize * reduce(operator.mul, data.shape[1:], 1)
    return min(data.shape[0], max(1, (2 ** 20) // rowsize))




def calculate_ld(nt_map_file, kgenomes_file, output_folder, window_size):
    """
    Calculate LD in windows for a reference genome dataset for a given set of SNPIds that are defined in the genotype_file
    """
    log.info('Calculating LD')
    # Load 1K genome
    kg_h5f = h5py.File(kgenomes_file, 'r')

    # load map file.
    with open(nt_map_file, 'rb') as f:
        snp_map_dict = pickle.load(f, encoding='latin1')

    # Figure out overlap (all genotype SNPs should be in the 1K genomes data)..
    for chrom in range(1, 23):
        log.info('Working on Chromosome %s' % chrom)
        chrom_str1 = 'chr%s' % chrom
        kg_cg = kg_h5f[chrom_str1]
        kg_sids = kg_cg['snp_ids'][...]
        chrom_dict = snp_map_dict[chrom_str1]
        g_sids = chrom_dict['sids']

        kg_filter = sp.in1d(kg_sids, g_sids)

        assert sp.sum(kg_filter) == len(g_sids), '..bug...'
        assert sp.all(kg_sids[kg_filter] == g_sids), '...bug'

        snps = kg_cg['snps'][...]
        snps = snps.compress(kg_filter, axis=0)

        snp_stds = kg_cg['snp_stds'][...]
        snp_stds = snp_stds.compress(kg_filter, axis=0)

        snp_means = kg_cg['snp_means'][...]
        snp_means = snp_means.compress(kg_filter, axis=0)

        norm_snps = sp.array((snps - snp_means) / snp_stds, dtype='single')

        # Iterate over SNPs and calculate LD
        num_snps, num_indivs = snps.shape

        ld_mats = []
        boundaries = []

        for snp_i in range(num_snps):
            start_i = max(0, snp_i - window_size / 2)
            end_i = min(snp_i + (window_size / 2) + 1, num_snps)

            X = norm_snps[start_i:end_i]
            D = sp.dot(X, X.T) / num_indivs

            ld_mats.append(D)
            boundaries.append([start_i, end_i])

        ld_dict = {'Ds':ld_mats, 'boundaries':boundaries, 'snp_means':snp_means, 'snp_stds':snp_stds, 'window_size':window_size}
        # Store things

        ld_file = '%s/LD' % output_folder + '_' + chrom_str1 + '.pickled.gz'
        log.info('Saving LD in %s' % ld_file)
        with gzip.open(ld_file, 'w') as f:
            pickle.dump(ld_dict, f, protocol=2)


def impute(genotype_file, ld_folder, output_file, validation_missing_rate=0.02, min_ld_r2_thres=0.02,
           regularization_factor=0.02, **kwargs):
    """
    Impute the missing SNPs from the reference genome dataset into a given genotype

    validation_missing_rate: The fraction of SNPs used to estimate the imputation accuracy.  2% seems enough to get SE<1%.
                             A smaller number will increase speed.

    min_ld_r2_thres: A minimum LD r2 value for SNPs to be used to impute from.  SNPs with large R2 are more informative for the imputation.
                     SNPs with r2 values close to 0 are effectively inconsequential for the imputation and can be left out
                     (which also speeds up the imputation).  Default is 0.02.
                     
    regularization_factor: It is a number bewtween 0 and 1, but should probably generally be small 0-0.1.  It effectively accounts 
                           for genotype and LD estimate errors, and increases the likelyhood that the LD matrix is invertible.
    """
    
    log.info('Starting imputation for %s using a missing rate of %s and minimum ld threshold of %s' % (genotype_file, validation_missing_rate, min_ld_r2_thres));
    g_h5f = h5py.File(genotype_file, 'r')
    imputed_snps_dict = {}

    pred_snps = []
    true_snps = []
    result = {'chr_stats':{}}
    
    log_extra = kwargs.get('log_extra', {'progress':0})
    if 'max_progress' not in log_extra:
        log_extra['max_progress'] = 100
    partial_progress_inc = (log_extra['max_progress'] - log_extra['progress']) / 22
    
    for chrom in range(1, 23):
        log_extra['progress'] += partial_progress_inc
        log.info('Working on Chromosome %d' % chrom, extra=log_extra)

        # Loading pre-calculated LD matrices (Note that these could perhaps be stored more efficiently.)
        chrom_str = 'Chr%d' % chrom
        with gzip.open('%s/LD' % ld_folder + '_' + chrom_str.lower() + '.pickled.gz', 'rb') as f:
            if version_info[0] < 3:
                ld_dict = pickle.load(f)
            else:
                ld_dict = pickle.load(f, encoding='latin1')

        g_cg = g_h5f[chrom_str]

        # Loading data
        snps = g_cg['snps'][...]
        Ds = ld_dict['Ds']
        snp_means = ld_dict['snp_means']
        snp_stds = ld_dict['snp_stds']

        # The snp vector to be returned
        imputed_snps = snps.copy()

        num_snps = len(snps)
        assert len(Ds) == num_snps, '..bug'
        num_snps_imputed = 0
        for snp_i in range(num_snps):
            
            if random.random() < validation_missing_rate and snps[snp_i] != -9:
                # Picking random SNPs with genotype information to estimate imputation accuracy.
                true_snp = snps[snp_i]
                snps[snp_i] = -9
            else:
                true_snp = -9
                
            if snps[snp_i] == -9:
                # Pull out LD matrix
                D = Ds[snp_i]

                # Determining the boundaries of the region.
                boundaries = ld_dict['boundaries'][snp_i]
                # start_i = max(0,snp_i-window_size/2)
                # end_i = min(snp_i+(window_size/2)+1,num_snps)
                start_i = boundaries[0]
                end_i = boundaries[1]

                # Obtaining the SNPs in the region, on which the imputation (together with LD) is based.
                reg_snps = snps[start_i:end_i]
                reg_snp_means = snp_means[start_i:end_i]
                reg_snp_means = reg_snp_means.flatten()
                reg_snp_stds = snp_stds[start_i:end_i]
                reg_snp_stds = reg_snp_stds.flatten()

                # The LD vector for the SNP to be imputed.
                loc_i = snp_i - start_i
                D_i = D[loc_i]

                # Filter SNPs that have missing genotypes
                ok_filter = reg_snps != -9

                # Filtering SNPs that are not in LD with the SNP to be imputed.  This saves time and may improve accuracy.
                ld_filter = (D_i ** 2 > min_ld_r2_thres)
                if sp.any(ok_filter * ld_filter):
                    ok_filter = ok_filter * ld_filter

                assert sp.sum(ok_filter) < len(reg_snps), '..bug'


                # Filtering the LD matrix.
                ok_D = (D[ok_filter])[:, ok_filter]
                ok_D_i = D_i[ok_filter]

                # Impute genotype.
                ok_D_inv = linalg.pinv((1 - regularization_factor) * ok_D + regularization_factor * sp.eye(len(ok_D)))  
                if sp.any(sp.isnan(ok_D_inv)):
                    log.warn('Matrix inversion failed!!')
                    log.warn('Setting SNP genotype to 1')
                    imputed_snp = 1
                else:
                    # Filtering the genotypes in the LD region.
                    ok_reg_snps = reg_snps[ok_filter]
                    ok_reg_snp_means = reg_snp_means[ok_filter]
                    ok_reg_snp_stds = reg_snp_stds[ok_filter]
                    ok_reg_norm_snps = (ok_reg_snps - ok_reg_snp_means) / ok_reg_snp_stds
#                     imputed_snp = sp.dot(ok_D_i,sp.dot(ok_D_inv,ok_reg_snps))  #A bug?
                    imputed_snp = sp.dot(ok_D_i, sp.dot(ok_D_inv, ok_reg_norm_snps))

                    # Transform imputed genotype to 0-2 scale
                    snp_mean = snp_means[snp_i][0]
                    snp_std = snp_stds[snp_i][0]
                    imputed_snp = imputed_snp * snp_std + snp_mean
                    if imputed_snp < 0:
                        imputed_snp = 0
                    elif imputed_snp > 2:
                        imputed_snp = 2
                
                if true_snp != -9:
                    # Estimate prediction accuracy
                    pred_snps.append(imputed_snp)
                    true_snps.append(true_snp)
                else:
                    # Counting the imputed SNPs with actual missing genotype information and setting the imputed value
                    num_snps_imputed += 1
                    # Storing imputed genotypes
                    imputed_snps[snp_i] = imputed_snp
                    
        result['chr_stats'][chrom_str] = num_snps_imputed
        log.info('Number of SNPs imputed so far: %d' % num_snps_imputed)
        imputed_snps_dict[chrom_str] = imputed_snps

    pred_r2 = (sp.corrcoef(pred_snps, true_snps)[0, 1]) ** 2
    log.info('Estimated prediction accuracy (R2): %0.4f' % pred_r2)
    result['pred_r2'] = float(pred_r2)

    if output_file:
        log.info('Writing imputed genotypes to file %s' % output_file)
        oh5f = h5py.File(output_file, 'w')
        for chrom in range(1, 23):
            log.info('Working on Chromosome %d' % chrom)

            chrom_str = 'Chr%d' % chrom
            g_cg = g_h5f[chrom_str]
            imputed_snps = imputed_snps_dict[chrom_str]

            # Loading data 
            sids = g_cg['sids'][...]
            nts = g_cg['nts'][...]
            positions = g_cg['positions'][...]

            cg = oh5f.create_group(chrom_str)
            cg.create_dataset('snps', data=imputed_snps, compression='lzf', chunks=True)
            cg.create_dataset('sids', data=sids, compression='lzf', chunks=True)
            cg.create_dataset('positions', data=positions, compression='lzf', chunks=True)
            cg.create_dataset('nts', data=nts, compression='lzf', chunks=True)
        oh5f.attrs['source'] = g_h5f.attrs['source']
        oh5f.attrs['version'] = g_h5f.attrs['version']
        oh5f.attrs['gender'] = g_h5f.attrs['gender']
        oh5f.close()
    g_h5f.close()
    return result



# def _window_size_plot_():
#     """
#     For debugging purposes
#     """
#     import pylab
#     pred_r2s = []
#     window_sizes = [4, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
#     for window_size in window_sizes:
# #         calc_ld(repos_dir+'imputor/tests/data/test_out_genotype.hdf5', repos_dir+'imputor/tests/data/ld_dict',window_size=window_size)
#         d = impute_23_and_genome(window_size=window_size)
#         pred_r2s.append(d['pred_r2'])
#     print pred_r2s
#     print window_sizes
#      
#     pylab.plot(window_sizes, pred_r2s, alpha=0.6)
#     pylab.ylabel('Prediction accuracy (R2)')
#     pylab.xlabel('Imputation LD window-size')
#     pylab.savefig(cloud_dir + 'tmp/tmp.png')


    
# For debugging purposes
# if __name__=='__main__':
#     Filter related indivs
#     gen_unrelated_eur_1k_data()
#     
#     prepare_nt_coding_key(cloud_dir+'Data/1Kgenomes/1K_genomes_v3_EUR_unrelated2.hdf5',
#                           repos_dir+'imputor/tests/data/test_genotype.hdf5',
#                           cloud_dir+'tmp/nt_map.pickled')
#     parse_hdf5_genotype(repos_dir+'imputor/tests/data/test_genotype.hdf5',
#                          cloud_dir+'tmp/nt_map.pickled',
#                          repos_dir+'imputor/tests/data/test_out_genotype.hdf5')
#     
#     window_size = int(sys.argv[1])
#     
#     calc_ld(cloud_dir+'tmp/nt_map.pickled', repos_dir+'imputor/tests/data/ld_dict',window_size=window_size)
#     impute_23_and_genome(genotype_file=cloud_dir+'tmp/1k_ind_4.hdf5',window_size=window_size)
#      window_size_plot()




