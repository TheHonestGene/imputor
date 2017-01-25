"""
Code to generate the nucleotide map for 23andme/Ancestry genotypes, for a subset of SNPs (which depends on the genotype data version).
"""
import logging
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

# MAYBE NOT NEEDED if it is done in prepare_haptmap_for_ld_calculate function
def create_coding_key_map(K_genomes_file, genotype_file, nt_map_file):
    """
    Creates a nucleotides map for the reference SNPs

    """
    log.info('Generating NT map')
    gf = h5py.File(genotype_file, 'r')
    kgf = h5py.File(K_genomes_file, 'r')
    snp_map_dict = {}
    num_snps = 0
    for chrom in range(1, 23):
        log.info('Working on chromosome %d' % chrom)
        kg_chrom_str = 'chr%d' % chrom
        chrom_str = 'Chr%d' % chrom

        # Get SNPs from genotype
        cg = gf[chrom_str]
        sids = cg['ids'][...]

        # Get SNP IDs from 1K genomes
        kcg = kgf[kg_chrom_str]
        kg_sids = kcg['snp_ids'][...]

        # Determine overlap between SNPs..
        kg_filter = sp.in1d(kg_sids, sids)
        kg_sids = kg_sids.compress(kg_filter, axis=0)
        kg_nts = kcg['nts'][...].compress(kg_filter, axis=0)
        kg_positions = kcg['positions'][...].compress(kg_filter, axis=0)

        # Check that nt are ok in genotype data, otherwise filter.
        sid_nt_map = {}
        positions = []
        ok_sids = []
        nts = []
        snp_i = 0
        for sid, kg_nt, kg_pos in zip(kg_sids, kg_nts, kg_positions):
            if tuple(kg_nt) not in ambig_nts:
                # All possible (and allowed) nucleotides strings
                ntm = {}
                a = kg_nt[0].decode()
                b = kg_nt[1].decode()
                ntm['--'] = -9
                ntm['-' + a] = -9
                ntm['-' + b] = -9
                ntm[a + '-'] = -9
                ntm[b + '-'] = -9
                ntm[a + a] = 0
                ntm[b + a] = 1
                ntm[a + b] = 1
                ntm[b + b] = 2
                if sid in sid_nt_map:
                    log.warn('Warning %s SNP is duplicate' % sid) 
                sid_nt_map[sid] = {'ntm':ntm, 'snp_i':snp_i}
                positions.append(kg_pos)
                nts.append(kg_nt)
                ok_sids.append(sid)
                snp_i += 1

        num_snps += len(sid_nt_map)

        # Sorting SNPs by positions
        sort_indices = sp.argsort(positions)
        if not sp.all(sort_indices == sp.arange(len(sort_indices))):
            positions = positions[sort_indices]
            ok_sids = ok_sids[sort_indices]
            nts = nts[sort_indices]

        snp_map_dict[kg_chrom_str] = {'sid_nt_map':sid_nt_map, 'positions':positions, 'nts':nts, 'sids':ok_sids}

    log.info('Found %d SNPs' % num_snps)
    log.info('Writing to file')
    with open(nt_map_file, 'wb') as f:
        pickle.dump(snp_map_dict, f, protocol=2)


