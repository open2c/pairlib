import numpy as np
import pandas as pd

from ._regions import assign_regs_c

def geomprog(factor, start=1):
    yield start
    while True:
        start *= factor
        yield start

def _geomrange(start, end, factor, endpoint):
    prev = np.nan
    for i in geomprog(factor, start):
        x = int(round(i))
        
        if x > end:
            break

        if x == prev:
            continue
        
        prev = x
        yield x

    if endpoint and prev != end:
        yield end

def geomrange(start, end, factor, endpoint=False):
    return np.fromiter(_geomrange(start, end, factor, endpoint), dtype=int)

def geomspace(start, end, num=50, endpoint=True):
    factor = (end / start) ** (1 / num)
    return geomrange(start, end, factor, endpoint=endpoint)


def _to_float(arr_or_scalar):
    if np.isscalar(arr_or_scalar):
        return float(arr_or_scalar)
    else:
        return np.asarray(arr_or_scalar).astype(float)

def assign_regs(chroms, pos, regs):
    gb_regs = regs.sort_values(['chrom', 'start', 'end']).groupby(['chrom'])
    
    regs_dict = {
        chrom.encode() : regs_per_chrom[['start','end']].values.flatten()
        for chrom, regs_per_chrom in gb_regs
    }
    
    return assign_regs_c(
        np.asarray(chroms).astype('bytes'),
        np.asarray(pos), 
        regs_dict)
    

def bins_pairs_by_distance(
    pairs_df, 
    dist_bins,
    regions=None,
    chromsizes=None,
    ):
   
    if regions is None:
        if chromsizes is None:
            region_starts1, region_starts2 = 0, 0
            region_ends1, region_ends2 = -1, -1
        else:
            region_starts1, region_starts2 = 0, 0
            region_ends1 = pairs_df.chrom1.map(chromsizes).fillna(1).astype(np.int64)
            region_ends2 = pairs_df.chrom2.map(chromsizes).fillna(1).astype(np.int64)
    else:
        _, region_starts1, region_ends1 = assign_regs(
            pairs_df.chrom1.values,
            pairs_df.pos1.values, 
            regions).T
        _, region_starts2, region_ends2 = assign_regs(
            pairs_df.chrom2.values,
            pairs_df.pos2.values, 
            regions).T        
    
     
    dist_bin_idxs = np.searchsorted(
        dist_bins, 
        pairs_df.eval('abs(pos1-pos2)'),
        side='right'
    )
    
    min_dist = np.where(dist_bin_idxs>0, 
                        dist_bins[dist_bin_idxs-1],
                        0)
    
    max_dist = np.where(dist_bin_idxs<len(dist_bins), 
                        dist_bins[dist_bin_idxs], 
                        np.iinfo(np.int64).max)
    
    is_same_region = (
        (pairs_df.chrom1.values == pairs_df.chrom2.values)
        & (region_starts1 == region_starts2)
    )
    min_dist = np.where(is_same_region, min_dist, 0)
    max_dist = np.where(is_same_region, max_dist, 0)
    
    pairs_reduced_df = pd.DataFrame(
        {'min_dist': min_dist,
         'max_dist': max_dist,
         'chrom1':pairs_df.chrom1.values,
         'chrom2':pairs_df.chrom2.values,
         'region_start1':region_starts1,
         'region_end1':region_ends1,
         'region_start2':region_starts2,
         'region_end2':region_ends2,
         'strand1':pairs_df.strand1,
         'strand2':pairs_df.strand2,
         'n_pairs':1
         },
        copy=False)

    pairs_reduced_gb = pairs_reduced_df.groupby(
        by=['chrom1','chrom2', 
            'region_start1', 'region_end1', 'region_start2', 'region_end2', 
            'strand1', 'strand2', 
            'min_dist', 'max_dist'])
    return pairs_reduced_gb.count()


def contact_areas_same_reg(
    min_dist, 
    max_dist, 
    region_length
    ):
    
    min_dist = _to_float(min_dist)
    max_dist = _to_float(max_dist)
    scaffold_length = _to_float(region_length)
    outer_areas = np.maximum(region_length - min_dist, 0) ** 2
    inner_areas = np.maximum(region_length - max_dist, 0) ** 2
    return 0.5 * (outer_areas - inner_areas) 

def _contact_areas_diff_reg(
    min_dist, 
    max_dist, 
    region_start1,
    region_end1,
    region_start2,
    region_end2
    ):
    
    return (contact_areas_same_reg(min_dist, max_dist, np.abs(region_end2 - region_start1))
            + contact_areas_same_reg(min_dist, max_dist, np.abs(region_end1 - region_start2))
            - contact_areas_same_reg(min_dist, max_dist, np.abs(region_start1 - region_start2))
            - contact_areas_same_reg(min_dist, max_dist, np.abs(region_end1 - region_end2))
           )
    
def _contact_areas_trans(
    min_dist,
    max_dist,
    region_length1,
    region_length2
    ):
    
    return (
        contact_areas_same_reg(min_dist, max_dist, region_length1+region_length2)
        -contact_areas_same_reg(min_dist, max_dist, region_length1)
        -contact_areas_same_reg(min_dist, max_dist, region_length2)
    )


def compute_scaling(
    pairs_df, 
    regions=None,
    chromsizes=None,
    dist_range=(int(1e1), int(1e9)), 
    n_dist_bins=8*8,
    ):

    dist_bins = geomspace(dist_range[0],dist_range[1],n_dist_bins)

    sc = bins_pairs_by_distance(
        pairs_df, 
        dist_bins,
        regions=regions,
        chromsizes=chromsizes
        )
    
    sc.reset_index(inplace=True)
    
#         if not (isinstance(regions, pd.DataFrame) and
#                  (set(regions.columns) == set(['chrom', 'start','end']))):
#             raise ValueError('regions must be provided as a dict or chrom-indexed Series of chromsizes or as a bedframe.')
            
    sc['n_bp2'] = contact_areas_same_reg(
        sc['min_dist'],
        sc['max_dist'],
        sc['region_end1'] - sc['region_start1']
        )
        
    return sc
