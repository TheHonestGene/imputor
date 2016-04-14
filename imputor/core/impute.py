"""
Code to impute the 23andme genome for the necessary SNPs.

"""
import h5py
import sys
import logging
import scipy as sp
from scipy import linalg
import gzip
import random
import io
import sn
try:
    from itertools import izip as zip
except ImportError: # will be 3.x series
    pass
import gzip
try:
    import cPickle as pickle
except ImportError: # will be 3.series
    import pickle
import numpy as np
import operator


ambig_nts = set([(b'A', b'T'), (b'T', b'A'), (b'G', b'C'), (b'C', b'G')])

log = logging.getLogger(__name__)

def _get_chunk_length(data):
    rowsize = data.dtype.itemsize * reduce(operator.mul,data.shape[1:],1)
    return min(data.shape[0],max(1,(2**20) // rowsize))


def convert_genotype_to_hdf5(csv_content,output_file,source=None):
    log.info('Convert genotype from text format to HDF5 format %s' % (output_file))
    # guess the source
    fd = io.StringIO(csv_content)
    if source is None:
        source = sn.guess_source_from_content(csv_content)
    version = ''
    snps = sn.parse(fd,source)
    start_chr = None
    f = h5py.File(output_file,'w')
    pos_snps =[]
    num_snps = 0
    for snp in snps:
        num_snps +=1
        if snp.chromosome != start_chr:
            if start_chr is not None:
                sorted(pos_snps, key=lambda x: x[0])
                positions,ids,snps = zip(*pos_snps)
                group.create_dataset('ids',(len(positions),),chunks=True,compression='lzf',dtype='S15',data=ids)
                group.create_dataset('positions',(len(positions),),chunks=True,compression='lzf',dtype='i8',data=positions)
                group.create_dataset('snps',(len(positions),),chunks=True,compression='lzf',dtype='S2',data=snps)
                pos_snps =[]
            start_chr = snp.chromosome
            group = f.create_group("Chr%s" % start_chr)
        pos_snps.append((int(snp.position),snp.name.encode('utf-8'),snp.genotype.encode('utf-8')))
    sorted(pos_snps, key=lambda x: x[0])
    positions,ids,snps = zip(*pos_snps)
    group.create_dataset('ids',(len(positions),),chunks=True,compression='lzf',dtype='S15',data=ids)
    group.create_dataset('positions',(len(positions),),chunks=True,compression='lzf',dtype='i8',data=positions)
    group.create_dataset('snps',(len(positions),),chunks=True,compression='lzf',dtype='S2',data=snps)
    # find out the version
    f.attrs['source'] = source
    # based on https://en.wikipedia.org/wiki/23andMe
    if source == '23andme':
        if num_snps <= 576000:
            version = 'v1'
        elif num_snps <= 597000:
            version = 'v2'
        # TODO current workaround because oauth genotype are > 1000000
        elif num_snps <= 611000 or num_snps > 1000000:
            version = 'v4'
        elif num_snps <= 992000:
            version = 'v3'
    f.attrs['version'] = version
    f.close()

def prepare_hapmap_for_ld_calculation(input_file,output_file):
    """
    Removes non-europeans and related individuals and monomorphic or unknown SNPs
    """
    h5f = h5py.File(input_file,'r')
    eur_filter = h5f['indivs']['continent'][...]=='EUR'
    num_indivs = sp.sum(eur_filter)
    K = sp.zeros((num_indivs,num_indivs), dtype='single')
    num_snps = 0
    log.info('Calculating kinship')
    for chrom in range(1,23):
        if chrom == 23:
            chrom='X'
        log.info('Working on Chromosome %s'%chrom)
        chrom_str = 'chr%s'%chrom
        log.debug('Loading SNPs')
        snps = h5f[chrom_str]['calldata/snps'][...]
        #filter non-europeans.
        log.debug('Filtering non-European individuals')
        snps = snps.compress(eur_filter,axis=1)
        log.debug('Filtering monomorphic SNPs')
        snp_stds = sp.std(snps,1)
        mono_morph_filter = snp_stds>0
        snps = snps.compress(mono_morph_filter,axis=0)
        snp_stds = snp_stds.compress(mono_morph_filter)
        log.debug('Normalizing SNPs')
        snp_means = sp.mean(snps,1)
        norm_snps = (snps - snp_means[sp.newaxis].T)/snp_stds[sp.newaxis].T
        log.debug('Updating kinship')
        K += sp.dot(norm_snps.T,norm_snps)
        num_snps += len(norm_snps)
    K = K/float(num_snps)
    log.info('Kinship calculation done using %d SNPs'%num_snps)

    #Filter individuals
    log.info('Filtering individuals')
    keep_indiv_set = set(range(num_indivs))
    for i in range(num_indivs):
        if i in keep_indiv_set:
            for j in range(i+1,num_indivs):
                if K[i,j]>0.05:
                    if j in keep_indiv_set:
                        keep_indiv_set.remove(j)
    keep_indivs = list(keep_indiv_set)
    keep_indivs.sort()
    log.info('Retained %d individuals' %len(keep_indivs))

    #Store in new file
    log.info('Now storing data.')
    oh5f = h5py.File(output_file,'w')
    indiv_ids = h5f['indivs']['indiv_ids'][eur_filter]
    indiv_ids = indiv_ids[keep_indivs]
    oh5f.create_dataset('indiv_ids',data =indiv_ids)
    for chrom in range(1,23):
        if chrom == 23:
            chrom='X'
        log.info('Working on Chromosome %s'%chrom)
        chrom_str = 'chr%s'%chrom
        snps = h5f[chrom_str]['calldata/snps'][...]
        #filter non-europeans.
        snps = snps.compress(eur_filter,axis=1)
        #Filter related
        snps = snps.compress(keep_indivs,axis=1)
        #filter monomorphic SNPs
        snp_stds = sp.std(snps,1)
        mono_morph_filter = snp_stds>0
        snps = snps.compress(mono_morph_filter,axis=0)
        #filter SNPs w missing NT values
        length = len(h5f[chrom_str]['variants/REF'])
        nts= np.hstack((h5f[chrom_str]['variants/REF'][:].reshape(length,1),h5f[chrom_str]['variants/ALT'][:].reshape(length,1)))
        nts = nts.compress(mono_morph_filter,axis=0)

        cg = oh5f.create_group(chrom_str)
        snp_chunk_length = _get_chunk_length(snps)
        cg.create_dataset('snps',data=snps,compression='lzf',chunks=(snp_chunk_length,snps.shape[1]))

        snp_stds = snp_stds.compress(mono_morph_filter)
        snp_means = sp.mean(snps,1)
        snps_means_trans = snp_means[sp.newaxis].T
        snps_stds_trans = snp_stds[sp.newaxis].T
        cg.create_dataset('snp_means',data=snps_means_trans,compression='lzf',chunks=(_get_chunk_length(snps_means_trans),1))
        cg.create_dataset('snp_stds',data=snps_stds_trans,compression='lzf',chunks=(_get_chunk_length(snps_stds_trans),1))

        snp_ids = h5f[chrom_str]['variants/ID'][...]
        snp_ids = snp_ids.compress(mono_morph_filter)
        cg.create_dataset('snp_ids',data=snp_ids,compression='lzf',chunks=(_get_chunk_length(snp_ids),))

        positions = h5f[chrom_str]['variants/POS'][...]
        positions = positions.compress(mono_morph_filter)
        cg.create_dataset('positions',data=positions,compression='lzf',chunks=(_get_chunk_length(positions),))
        cg.create_dataset('nts',data=nts,compression='lzf',chunks=(_get_chunk_length(nts),2))
    log.info('File successfully stored')
    oh5f.close()
    h5f.close()




# MAYBE NOT NEEDED if it is done in prepare_haptmap_for_ld_calculate function
def create_coding_key_map(K_genomes_file,genotype_file,nt_map_file):
    """
    Creates a nucleotides map for the reference SNPs

    """
    log.info('Generating NT map')
    gf = h5py.File(genotype_file,'r')
    kgf = h5py.File(K_genomes_file,'r')
    snp_map_dict = {}
    num_snps = 0
    for chrom in range(1,23):
        log.info('Working on chromosome %d' % chrom)
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
        kg_sids = kg_sids.compress(kg_filter,axis=0)
        kg_nts = kcg['nts'][...].compress(kg_filter,axis=0)
        kg_positions = kcg['positions'][...].compress(kg_filter,axis=0)

        #Check that nt are ok in genotype data, otherwise filter.
        sid_nt_map = {}
        positions = []
        ok_sids = []
        nts = []
        snp_i = 0
        for sid, kg_nt, kg_pos in zip(kg_sids, kg_nts, kg_positions):
            snp = sid_dict[sid]
            if tuple(kg_nt) not in ambig_nts:
                # All possible (and allowed) nucleotides strings
                ntm = {}
                a = kg_nt[0].decode()
                b = kg_nt[1].decode()
                ntm['--']=-9
                ntm['-'+a]=-9
                ntm['-'+b]=-9
                ntm[a+'-']=-9
                ntm[b+'-']=-9
                ntm[a+a]=0
                ntm[b+a]=1
                ntm[a+b]=1
                ntm[b+b]=2
                sid_nt_map[sid]={'ntm':ntm, 'snp_i':snp_i}
                positions.append(kg_pos)
                nts.append(kg_nt)
                ok_sids.append(sid)
                snp_i += 1
        num_snps += len(sid_nt_map)

        #Sorting SNPs by positions
        sort_indices = sp.argsort(positions)
        if not sp.all(sort_indices==sp.arange(len(sort_indices))):
            positions = positions[sort_indices]
            ok_sids = ok_sids[sort_indices]
            nts = nts[sort_indices]

        snp_map_dict[kg_chrom_str]={'sid_nt_map':sid_nt_map, 'positions':positions, 'nts':nts, 'sids':ok_sids}

    log.info('Found %d SNPs'%num_snps)
    log.info('Writing to file')
    with open(nt_map_file, 'wb') as f:
        pickle.dump(snp_map_dict, f, protocol=2)


def convert_genotype_nt_key_encoding(input_file,output_file,nt_map_file,**kwargs):
    """
    Convert the SNPs from nt form to numeric form in genotype to be imputed
    """
    log_extra = kwargs.get('log_extra',{'progress':0})
    if 'max_progress' not in log_extra:
        log_extra['max_progress'] = 100
    partial_progress_inc = (log_extra['max_progress']-log_extra['progress'])/22
    
    
    log.info('Loading NT map from file: %s'%nt_map_file)
    with open(nt_map_file, 'rb') as f:
        snp_map_dict = pickle.load(f,encoding='latin1')

    log.info('Parsing individual genotype: %s'%input_file)
    h5f = h5py.File(input_file,'r')
    #prepare output file
    oh5f = h5py.File(output_file,'w')
    try:
        tot_num_parsed_snps = 0
        result = {'total_num_parsed_snps':tot_num_parsed_snps,'chr_stats':{}}
        for chrom in range(1,23):
            log_extra['progress']+=partial_progress_inc
            log.info('Working on chromosome %s'%chrom,extra=log_extra)
            kg_chrom_str = 'chr%s'%chrom
            chrom_str = 'Chr%s'%chrom
            cg = h5f[chrom_str]
            sids = cg['ids'][...]
            raw_snps = cg['snps'][...]

            #Get the nucleotides coding map (from 1K genomes project).
            chrom_dict = snp_map_dict[kg_chrom_str]
            sid_nt_map = chrom_dict['sid_nt_map']
            n = len(sid_nt_map)
            snps = sp.repeat(-9, n) #Creating the SNP with fixed size
            num_not_found = 0
            num_misunderstood = 0
            num_parsed_ok = 0
            for sid, nt in zip(sids,raw_snps):
                try:
                    d = sid_nt_map[sid]
                except Exception:
                    num_not_found +=1
                    continue
                try:
                    nt_val = d['ntm'][nt.decode()] # decode required for py3
                except Exception:
                    num_misunderstood +=1
                    continue
                snps[d['snp_i']] = nt_val
                num_parsed_ok += 1


            log.info("%d SNPs weren't found and %d SNPs had unrecognizable nucleotides"%(num_not_found,num_misunderstood))
            log.info("%d SNPs were parsed ok."%num_parsed_ok)
            tot_num_parsed_snps +=num_parsed_ok
            result['chr_stats'][chrom_str] = {'num_not_found':num_not_found,'num_misunderstood':num_misunderstood,'num_parsed_ok':num_parsed_ok}
            #Not sure what information we need, perhaps only the SNPs?

            assert len(snps)==len(chrom_dict['sids'])==len(chrom_dict['positions'])==len(chrom_dict['nts']), '..bug'
            #Store information
            sids = chrom_dict['sids']
            positions = chrom_dict['positions']
            nts = chrom_dict['nts']
            cg = oh5f.create_group(chrom_str)
            cg.create_dataset('snps', data=snps,compression='lzf',chunks=True)
            cg.create_dataset('sids', data=chrom_dict['sids'],compression='lzf',chunks=True)
            cg.create_dataset('positions', data=chrom_dict['positions'],compression='lzf',chunks=True)
            cg.create_dataset('nts', data=chrom_dict['nts'],compression='lzf',chunks=True)
        log.info('In total %d SNPs were parsed.'%tot_num_parsed_snps)
        result['total_num_parsed_snps'] = tot_num_parsed_snps
        oh5f.attrs['source'] = h5f.attrs['source']
        oh5f.attrs['version'] = h5f.attrs['version']
    finally:
        h5f.close()
        oh5f.close()
    return result
    #return genome_dict
    
    

def gen_unrelated_eur_1k_data(out_file='Data/1Kgenomes/1K_genomes_v3_EUR_unrelated2.hdf5'):
    h5f = h5py.File(cloud_dir+'Data/1Kgenomes/1K_genomes_v3.hdf5')
    eur_filter = h5f['indivs']['continent'][...]=='EUR'
    num_indivs = sp.sum(eur_filter)
    K = sp.zeros((num_indivs,num_indivs), dtype='single')
    num_snps = 0
    log.info('Calculating kinship')
    for chrom in range(1,23):
        log.debug('Working on Chromosome %d'%chrom)
        chrom_str = 'chr%d'%chrom
        log.debug('Loading SNPs')
        snps = h5f[chrom_str]['raw_snps'][...]
        #filter non-europeans.
        log.debug('Filtering non-European individuals')
        snps = snps[:,eur_filter]
        log.debug('Filtering monomorphic SNPs')
        snp_stds = sp.std(snps,1)
        mono_morph_filter = snp_stds>0
        snps = snps[mono_morph_filter]
        snp_stds = snp_stds[mono_morph_filter]
        log.debug('Filter SNPs with missing NT information')
        nts = h5f[chrom_str]['nts'][...]
        nts = nts[mono_morph_filter]
        nt_filter = sp.all(nts>0,1)
        snps = snps[nt_filter]
        snp_stds = snp_stds[nt_filter]
        log.debug('Normalizing SNPs')
        snp_means = sp.mean(snps,1)
        norm_snps = (snps - snp_means[sp.newaxis].T)/snp_stds[sp.newaxis].T
        log.debug('Updating kinship')
        K += sp.dot(norm_snps.T,norm_snps)
        num_snps += len(norm_snps)
    K = K/float(num_snps)
    log.info('Kinship calculation done using %d SNPs\n'%num_snps)

    
    #Filter individuals
    log.info('Filtering individuals')
    keep_indiv_set = set(range(num_indivs))
    for i in range(num_indivs):
        if i in keep_indiv_set:
            for j in range(i+1,num_indivs):
                if K[i,j]>0.05:
                    if j in keep_indiv_set:
                        keep_indiv_set.remove(j)
    keep_indivs = list(keep_indiv_set)
    keep_indivs.sort()
    log.info('Retained %d individuals\n'%len(keep_indivs))

    #Store in new file
    log.info('Now storing data.')
    oh5f = h5py.File(cloud_dir+out_file,'w')
    indiv_ids = h5f['indivs']['indiv_ids'][eur_filter]
    indiv_ids = indiv_ids[keep_indivs]
    oh5f.create_dataset('indiv_ids',data =indiv_ids)
    for chrom in range(1,23):
        log.debug('Working on Chromosome %d'%chrom)
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
    log.info('File successfully stored')
    oh5f.close()
    h5f.close()


def calculate_ld(nt_map_file,kgenomes_file, output_folder, window_size):
    """
    Calculate LD in windows for a reference genome dataset for a given set of SNPIds that are defined in the genotype_file
    """
    log.info('Calculating LD')
    #Load 1K genome
    kg_h5f = h5py.File(kgenomes_file,'r')

    #load map file.
    with open(nt_map_file, 'rb') as f:
        snp_map_dict = pickle.load(f,encoding='latin1')

    #Figure out overlap (all genotype SNPs should be in the 1K genomes data)..
    for chrom in range(1,23):
        log.info('Working on Chromosome %s'%chrom)
        chrom_str1 = 'chr%s'%chrom
        kg_cg = kg_h5f[chrom_str1]
        kg_sids = kg_cg['snp_ids'][...]
        chrom_dict = snp_map_dict[chrom_str1]
        g_sids = chrom_dict['sids']

        kg_filter = sp.in1d(kg_sids,g_sids)

        assert sp.sum(kg_filter)==len(g_sids), '..bug...'
        assert sp.all(kg_sids[kg_filter]==g_sids), '...bug'

        snps = kg_cg['snps'][...]
        snps = snps.compress(kg_filter,axis=0)

        snp_stds = kg_cg['snp_stds'][...]
        snp_stds = snp_stds.compress(kg_filter,axis=0)

        snp_means = kg_cg['snp_means'][...]
        snp_means = snp_means.compress(kg_filter, axis=0)

        norm_snps = sp.array((snps - snp_means)/snp_stds,dtype='single')

        #Iterate over SNPs and calculate LD
        num_snps,num_indivs = snps.shape

        ld_mats = []
        boundaries = []

        for snp_i in range(num_snps):
            start_i = max(0,snp_i-window_size/2)
            end_i = min(snp_i+(window_size/2)+1,num_snps)

            X = norm_snps[start_i:end_i]
            D =  sp.dot(X,X.T)/num_indivs

            ld_mats.append(D)
            boundaries.append([start_i,end_i])

        ld_dict={'Ds':ld_mats,'boundaries':boundaries, 'snp_means':snp_means, 'snp_stds':snp_stds,'window_size':window_size}
        #Store things

        ld_file = '%s/LD' % output_folder +'_'+chrom_str1+'.pickled.gz'
        log.info('Saving LD in %s' % ld_file)
        with gzip.open(ld_file,'w') as f:
            pickle.dump(ld_dict, f, protocol=2)


def impute(genotype_file,ld_folder,output_file,validation_missing_rate=0.02, min_ld_r2_thres = 0.02,**kwargs):
    """
    Impute the missing SNPs from the reference genome dataset into a given genotype

    validation_missing_rate: The fraction of SNPs used to estimate the imputation accuracy.  2% seems enough to get SE<1%.
                             A smaller number will increase speed.

    min_ld_r2_thres: A minimum LD r2 value for SNPs to be used to impute from.  SNPs with large R2 are more informative for the imputation.
                     SNPs with r2 values close to 0 are effectively inconsequential for the imputation and can be left out
                     (which also speeds up the imputation).  Default is 0.02.
    """
    log.info('Starting imputation for %s using a missing rate of %s and minimum ld threshold of %s' % (genotype_file,validation_missing_rate,min_ld_r2_thres));
    g_h5f = h5py.File(genotype_file,'r')
    imputed_snps_dict = {}

    pred_snps = []
    true_snps = []
    result = {'chr_stats':{}}
    
    log_extra = kwargs.get('log_extra',{'progress':0})
    if 'max_progress' not in log_extra:
        log_extra['max_progress'] = 100
    partial_progress_inc = (log_extra['max_progress']-log_extra['progress'])/22
    
    for chrom in range(1,23):
        log_extra['progress']+=partial_progress_inc
        log.info('Working on Chromosome %d'%chrom,extra=log_extra)

        #Loading pre-calculated LD matrices (Note that these could perhaps be stored more efficiently.)
        chrom_str = 'Chr%d'%chrom
        with gzip.open('%s/LD'%ld_folder +'_'+chrom_str.lower()+'.pickled.gz','rb') as f:
            ld_dict = pickle.load(f,encoding='latin1')

        g_cg = g_h5f[chrom_str]

        #Loading data
        snps = g_cg['snps'][...]
        Ds = ld_dict['Ds']
        snp_means = ld_dict['snp_means']
        snp_stds = ld_dict['snp_stds']

        #The snp vector to be returned
        imputed_snps = snps.copy()

        num_snps = len(snps)
        assert len(Ds)==num_snps,'..bug'
        num_snps_imputed = 0
        for snp_i in range(num_snps):

            if random.random()<validation_missing_rate and snps[snp_i] !=-9:
                #Picking random SNPs with genotype information to estimate imputation accuracy.
                true_snp = snps[snp_i]
                snps[snp_i] = -9
            else:
                true_snp=-9

            if snps[snp_i] ==-9:
                #Pull out LD matrix
                D = Ds[snp_i]

                #Determining the boundaries of the region.
                boundaries = ld_dict['boundaries'][snp_i]
                #start_i = max(0,snp_i-window_size/2)
                #end_i = min(snp_i+(window_size/2)+1,num_snps)
                start_i = boundaries[0]
                end_i = boundaries[1]

                #Obtaining the SNPs in the region, on which the imputation (together with LD) is based.
                reg_snps = snps[start_i:end_i]
                reg_snp_means = snp_means[start_i:end_i]
                reg_snp_means = reg_snp_means.flatten()
                reg_snp_stds = snp_stds[start_i:end_i]
                reg_snp_stds = reg_snp_stds.flatten()

                #The LD vector for the SNP to be imputed.
                loc_i = snp_i-start_i
                D_i = D[loc_i]

                #Filter SNPs that have missing genotypes
                ok_filter = reg_snps !=-9

                #Filtering SNPs that are not in LD with the SNP to be imputed.  This saves time and may improve accuracy.
                ld_filter = (D_i**2>min_ld_r2_thres)
                if sp.any(ok_filter*ld_filter):
                    ok_filter = ok_filter*ld_filter

                #Filtering the genotypes in the LD region.
                assert sp.sum(ok_filter)<len(reg_snps), '..bug'
                ok_reg_snps = reg_snps[ok_filter]
                ok_reg_snp_means = reg_snp_means[ok_filter]
                ok_reg_snp_stds = reg_snp_stds[ok_filter]
                ok_reg_norm_snps = (ok_reg_snps-ok_reg_snp_means) /ok_reg_snp_stds

                #Filtering the LD matrix.
                ok_D = (D[ok_filter])[:,ok_filter]
                ok_D_i = D_i[ok_filter]

                #Impute genotype.
                ok_D_inv = linalg.pinv(0.95*ok_D+0.05*sp.eye(len(ok_D)))
                if sp.any(sp.isnan(ok_D_inv)):
                    log.warn('Matrix inversion failed!!')
                    log.warn('Setting SNP genotype to 1')
                    imputed_snp=1
                else:
                    imputed_snp = sp.dot(ok_D_i,sp.dot(ok_D_inv,ok_reg_snps))


                    #Transform imputed genotype to 0-2 scale
                    snp_mean = snp_means[snp_i][0]
                    snp_std = snp_stds[snp_i][0]
                    imputed_snp = imputed_snp*snp_std+snp_mean
                    if imputed_snp<0:
                        imputed_snp=0
                    elif imputed_snp>2:
                        imputed_snp=2

                #Storing imputed genotypes
                imputed_snps[snp_i] = imputed_snp
                if true_snp!=-9:
                    #Estimate prediction accuracy
                    pred_snps.append(imputed_snp)
                    true_snps.append(true_snp)
                else:
                    #Counting the imputed SNPs with actual missing genotype information
                    num_snps_imputed += 1
        result['chr_stats'][chrom_str] = num_snps_imputed
        log.info('Number of SNPs imputed so far: %d'%num_snps_imputed)
        imputed_snps_dict[chrom_str] = imputed_snps

    pred_r2 = (sp.corrcoef(pred_snps, true_snps)[0,1])**2
    log.info('Estimated prediction accuracy (R2): %0.4f'%pred_r2)
    result['pred_r2'] = float(pred_r2)

    if output_file:
        log.info('Writing imputed genotypes to file %s' % output_file)
        oh5f = h5py.File(output_file,'w')
        for chrom in range(1,23):
            log.info('Working on Chromosome %d'%chrom)

            chrom_str = 'Chr%d'%chrom
            g_cg = g_h5f[chrom_str]
            imputed_snps = imputed_snps_dict[chrom_str]

            #Loading data 
            sids = g_cg['sids'][...]
            nts = g_cg['nts'][...]
            positions = g_cg['positions'][...]

            cg = oh5f.create_group(chrom_str)
            cg.create_dataset('snps', data=imputed_snps,compression='lzf',chunks=True)
            cg.create_dataset('sids', data=sids,compression='lzf',chunks=True)
            cg.create_dataset('positions', data=positions,compression='lzf',chunks=True)
            cg.create_dataset('nts', data=nts,compression='lzf',chunks=True)
        oh5f.attrs['source'] = g_h5f.attrs['source']
        oh5f.attrs['version'] = g_h5f.attrs['version']
        oh5f.close()
    g_h5f.close()
    return result
