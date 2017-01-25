"""
Code to parse 23andme/Ancestry.com genomes and store in an internal HDF5 format.

"""
from sys import version_info
import logging
import h5py
import scipy as sp
import sn
import io

try:
    from itertools import izip as zip
except ImportError:  # will be 3.x series
    pass
try:
    import cPickle as pickle
except ImportError:  # will be 3.series
    import pickle

log = logging.getLogger(__name__)


def convert_genotype_nt_key_encoding(input_file, output_file, nt_map_file, **kwargs):
    """
    Convert the SNPs from nt form to numeric form in genotype to be imputed
    """
    log_extra = kwargs.get('log_extra', {'progress':0})
    if 'max_progress' not in log_extra:
        log_extra['max_progress'] = 100
    partial_progress_inc = (log_extra['max_progress'] - log_extra['progress']) / 22
    
    
    log.info('Loading NT map from file: %s' % nt_map_file)
    with open(nt_map_file, 'rb') as f:
        if version_info[0] < 3:
            snp_map_dict = pickle.load(f)
        else:
            snp_map_dict = pickle.load(f, encoding='latin1')
            
    log.info('Parsing individual genotype: %s' % input_file)
    h5f = h5py.File(input_file, 'r')
    # prepare output file
    oh5f = h5py.File(output_file, 'w')
    try:
        tot_num_parsed_snps = 0
        result = {'total_num_parsed_snps':tot_num_parsed_snps, 'chr_stats':{}}
        for chrom in range(1, 23):
            log_extra['progress'] += partial_progress_inc
            log.info('Working on chromosome %s' % chrom, extra=log_extra)
            kg_chrom_str = 'chr%s' % chrom
            chrom_str = 'Chr%s' % chrom
            cg = h5f[chrom_str]
            sids = cg['ids'][...]
            raw_snps = cg['snps'][...]

            # Get the nucleotides coding map (from 1K genomes project).
            chrom_dict = snp_map_dict[kg_chrom_str]
            sid_nt_map = chrom_dict['sid_nt_map']
            n = len(sid_nt_map)
            snps = sp.repeat(-9, n)  # Creating the SNP with fixed size
            num_not_found = 0
            num_misunderstood = 0
            num_parsed_ok = 0
            for sid, nt in zip(sids, raw_snps):
                try:
                    d = sid_nt_map[sid]
                except Exception:
                    num_not_found += 1
                    continue
                try:
                    nt_val = d['ntm'][nt.decode()]  # decode required for py3
                except Exception:
                    num_misunderstood += 1
                    continue
                snps[d['snp_i']] = nt_val
                num_parsed_ok += 1
            log.info("%d SNPs weren't found and %d SNPs had unrecognizable nucleotides" % (num_not_found, num_misunderstood))
            log.info("%d SNPs were parsed ok." % num_parsed_ok)
            tot_num_parsed_snps += num_parsed_ok
            result['chr_stats'][chrom_str] = {'num_not_found':num_not_found, 'num_misunderstood':num_misunderstood, 'num_parsed_ok':num_parsed_ok}
            # Not sure what information we need, perhaps only the SNPs?

            assert len(snps) == len(chrom_dict['sids']) == len(chrom_dict['positions']) == len(chrom_dict['nts']), '..bug'
            # Store information
            cg = oh5f.create_group(chrom_str)
            cg.create_dataset('snps', data=snps, compression='lzf', chunks=True)
            cg.create_dataset('sids', data=chrom_dict['sids'], compression='lzf', chunks=True)
            cg.create_dataset('positions', data=chrom_dict['positions'], compression='lzf', chunks=True)
            cg.create_dataset('nts', data=chrom_dict['nts'], compression='lzf', chunks=True)
        log.info('In total %d SNPs were parsed.' % tot_num_parsed_snps)
        result['total_num_parsed_snps'] = tot_num_parsed_snps
        oh5f.attrs['source'] = h5f.attrs['source']
        oh5f.attrs['version'] = h5f.attrs['version']
        oh5f.attrs['gender'] = h5f.attrs['gender']
    finally:
        h5f.close()
        oh5f.close()
    return result
    # return genome_dict

def convert_genotype_to_hdf5(csv_content, output_file, source=None):
    log.info('Convert genotype from text format to HDF5 format %s' % (output_file))
    # guess the source
    fd = io.StringIO(csv_content)
    if source is None:
        source = sn.guess_source_from_content(csv_content)
    version = ''
    snps = sn.parse(fd, source)
    start_chr = None
    f = h5py.File(output_file, 'w')
    pos_snps = []
    num_snps = 0
    group = None  # Eliminates syntax error highlight in Eclipse.
    for snp in snps:
        num_snps += 1
        if snp.chromosome != start_chr:
            if start_chr is not None:
                sorted(pos_snps, key=lambda x: x[0])
                positions, ids, snps = zip(*pos_snps)
                group.create_dataset('ids', (len(positions),), chunks=True, compression='lzf', dtype='S15', data=ids)
                group.create_dataset('positions', (len(positions),), chunks=True, compression='lzf', dtype='i8', data=positions)
                group.create_dataset('snps', (len(positions),), chunks=True, compression='lzf', dtype='S2', data=snps)
                pos_snps = []
            start_chr = snp.chromosome
            group = f.create_group("Chr%s" % start_chr)
        pos_snps.append((int(snp.position), snp.name.encode('utf-8'), snp.genotype.encode('utf-8')))
    sorted(pos_snps, key=lambda x: x[0])
    positions, ids, snps = zip(*pos_snps)
    group.create_dataset('ids', (len(positions),), chunks=True, compression='lzf', dtype='S15', data=ids)
    group.create_dataset('positions', (len(positions),), chunks=True, compression='lzf', dtype='i8', data=positions)
    group.create_dataset('snps', (len(positions),), chunks=True, compression='lzf', dtype='S2', data=snps)
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
    f.attrs['gender'] = 'm' if 'ChrY' in f.keys() else 'f'
    f.close()