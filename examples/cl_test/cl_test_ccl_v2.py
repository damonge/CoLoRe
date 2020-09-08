'''
    Module built to compute the data needed to make the CLs analysis comparing with theoretical expectations. 

    It can be imported as a module or used as script.
'''

import logging
log = logging.getLogger(__name__)
from contextlib import contextmanager
import warnings

import numpy as np
import healpy as hp
import pyccl as ccl
import os, sys
from astropy.io import fits
from itertools import combinations_with_replacement
from datetime import datetime

import argparse
import json

def savetofile(location,variables,variables_names):
    for i,variable in enumerate(variables):
        np.savetxt( location + '/' + variables_names[i] + '.dat', variable)
    return
    
@contextmanager   # Code to avoid output temporarily
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout
        
def fxn():
    warnings.warn("deprecrated", DeprecationWarning)
    
    
def compute_all_cls(sim_path, source=1, nside=128, max_files=None, downsampling=1, zbins=[-1,0.15,1], nz_h = 50, nz_min=0, nz_max=None):
    '''Method to compute all cls from the anafast function using output from CoLoRe.

    Args:
        sim_path (str): Path where the CoLoRe simulation is located.
        source (int, optional): Source from which to compute data (default: 1)
        nside (int, optional): nside to use (default:128)
        max_files (int, optional): number of srcs files to consider (default: None, consider all the files) 
        downsampling (float, optional): downsampling to apply to the data (from 0 to 1) (default: 1)
        zbins (array of floats, optional): defines the binning in redshift of the data analysis (default: [0,0.15,0.5])
        nz_h (int, optional): pixelization of the redshift analysis (default: 50)
        nz_min (float, optional): min redshift for the N(z) histogram (default: 0)
        nz_max (float, optional): max redshift for the N(z) histogram (default: None (set to the max value in zbins))
    
    Returns: 
        Tuple given by (shotnoise, pairs, nz_tot, z_nz, d_values, cl_dd_t, cl_dm_t, cl_mm_t):
            d_values: Array of shape (3,6,384) with all the pairs in [(0,0),(0,1),(1,1)] for all the combinations TT, TE, TB, EE, EB, BB for the anafast function. 
    '''

    # Hubble constant
    h = 0.7

    nside = nside
    npix = hp.nside2npix(nside)

    nbins   = len(zbins) - 1
    pairs   = list(combinations_with_replacement(range(nbins), r=2))


    # This will contain the N(z) of the different bins
    nz_tot = np.zeros([nbins, nz_h])
    nz_max = nz_max if nz_max is not None else zbins[-1]
    
    sigz = 0.03

    # These will be the density and ellipticity maps
    nmap = np.zeros([nbins, npix])
    dmap = np.zeros([nbins, npix])
    e1map = np.zeros([nbins, npix])
    e2map = np.zeros([nbins, npix])

    # Now we loop over all source files
    if max_files is None:
        file_limit = False
    else:
        file_limit = True

    ifile = 0

    log.info('Reading output files...\n')
    while os.path.isfile(sim_path + '/out_srcs_s1_%d.fits' % ifile) and ( (not file_limit) or ifile <= max_files):
        hdulist = fits.open(sim_path + '/out_srcs_s1_%d.fits' % ifile)
        d = hdulist[1].data
        n_g = len(d)

        # Generate random photo-z
        z_photo = d['Z_COSMO'] + sigz*(1+d['Z_COSMO'])*np.random.randn(n_g)

        if downsampling != 1: 
            d_mask = np.random.random(len(z_photo)) < downsampling #pylint: disable=no-member
        else:
            d_mask = np.full(len(z_photo), True, dtype=bool)

        masks = []
        for i in range(nbins):
            masks.append( (zbins[i] <= z_photo) & (z_photo < zbins[i+1]) )
        
        masks = [(mask & d_mask) for mask in masks]  # applying downsampling here

        # For each bin, add to the number and ellipticity maps
        for ibin, msk in enumerate( masks ):
            dd = d[msk]

            pix = hp.ang2pix(nside,
                            np.radians(90-dd['DEC']),
                            np.radians(dd['RA']))
            n = np.bincount(pix, minlength=npix)
            e1 = np.bincount(pix, minlength=npix, weights=dd['E1'])
            e2 = np.bincount(pix, minlength=npix, weights=dd['E2'])
            nmap[ibin, :] += n
            e1map[ibin, :] += e1
            e2map[ibin, :] += e2

            # Add also to N(z)
            nz, z_edges = np.histogram(dd['Z_COSMO'], bins=nz_h,
                                    range=[nz_min, nz_max])
            nz_tot[ibin, :] += nz
        ifile += 1

    # Midpoint of N(z) histogram
    z_nz = 0.5*(z_edges[1:] + z_edges[:-1])

    log.info('Computing maps...\n')
    # Compute <e> map and overdensity map
    # Compute also shot noise level
    shotnoise = np.zeros(nbins)
    for ib in range(nbins):
        ndens = (np.sum(nmap[ib])+0.0)/(4*np.pi)
        shotnoise[ib] = 1./ndens
        e1map[ib, :] = e1map[ib]/nmap[ib]
        e1map[ib, nmap[ib] <= 0] = 0
        e2map[ib, :] = e2map[ib]/nmap[ib]
        e2map[ib, nmap[ib] <= 0] = 0
        dmap[ib, :] = (nmap[ib, :] + 0.0) / np.mean(nmap[ib] + 0.0) - 1

    # Read P(k) theory prediction
    zs = []
    pks_dd = []
    pks_dm = []
    pks_mm = []
    z_pk = 0.
    log.info('Reading pk files...\n')
    while os.path.isfile(sim_path + "/out_pk_srcs_pop0_z%.3lf.txt" % z_pk):
        ks, pdd, pdm, pmm = np.loadtxt(sim_path + "/out_pk_srcs_pop0_z%.3lf.txt" % z_pk, unpack=True)
        # The delta-delta prediction involves some Fourier transforms that make it unstable
        # at high-k, so we just set it to zero.
        pdd[pmm<1E-30] = 0
        pks_dd.append(pdd)
        pks_dm.append(pdm)
        pks_mm.append(pmm)
        zs.append(z_pk)
        z_pk += 0.1
    # Reverse order (because CCL needs increasing scale factor - not redshift).
    # Also, CCL uses non-h units.
    zs = np.array(zs)[::-1]
    ks = np.array(ks) * h
    pks_dd = np.array(pks_dd)[::-1, :] / h**3
    pks_dm = np.array(pks_dm)[::-1, :] / h**3
    pks_mm = np.array(pks_mm)[::-1, :] / h**3

    log.info('Creating CCL structures...\n')
    # Create CCL P(k) structures
    cosmo = ccl.Cosmology(Omega_c=0.25, Omega_b=0.05, h=h, n_s=0.96, sigma8=0.8)
    pk2d_dd = ccl.Pk2D(a_arr=1./(1+zs), lk_arr=np.log(ks), pk_arr=pks_dd,
                    is_logp=False, extrap_order_hik=0, extrap_order_lok=0)
    pk2d_dm = ccl.Pk2D(a_arr=1./(1+zs), lk_arr=np.log(ks), pk_arr=pks_dm,
                    is_logp=False, extrap_order_hik=0, extrap_order_lok=0)
    pk2d_mm = ccl.Pk2D(a_arr=1./(1+zs), lk_arr=np.log(ks), pk_arr=pks_mm,
                    is_logp=False, extrap_order_hik=0, extrap_order_lok=0)
    # These can be evaluated as follows:
    pk2d_dd.eval(k=0.1, a=1/(1+0.2), cosmo=cosmo) # This was printed by David!

    # Create a number counts and a weak lensing tracer for each of the bins
    tr_d = [ccl.NumberCountsTracer(cosmo, False, (z_nz, nz_tot[i]), bias=(z_nz, np.ones_like(z_nz)))
            for i in range(nbins)]
    tr_l = [ccl.WeakLensingTracer(cosmo, (z_nz, nz_tot[i])) for i in range(nbins)]

    log.info('Computing CCL power spectra...\n')
    # Compute power spectra. I'm only doing delta-delta here.
    larr = np.arange(3*nside)
    cl_dd_t = np.array([ccl.angular_cl(cosmo, tr_d[p1], tr_d[p2], larr, p_of_k_a=pk2d_dd)
                        for p1, p2 in pairs])
    cl_dm_t = np.array([ccl.angular_cl(cosmo, tr_d[p1], tr_l[p2], larr, p_of_k_a=pk2d_dm)
                        for p1, p2 in pairs])
    cl_mm_t = np.array([ccl.angular_cl(cosmo, tr_l[p1], tr_l[p2], larr, p_of_k_a=pk2d_mm)
                        for p1, p2 in pairs])
    cl_md_t = np.array([ccl.angular_cl(cosmo, tr_d[p1], tr_l[p2], larr, p_of_k_a=pk2d_dm)
                        for p2, p1 in pairs])

    log.info('Computing values from data...\n')
    import pymaster as nmt

    b=nmt.NmtBin.from_edges(np.arange(3*nside+1)[:-1],
                            np.arange(3*nside+1)[1:])
    # Spin-0 fields
    mone = np.ones(npix)
    f0 = [nmt.NmtField(mone, [d], n_iter=0) for d in dmap]
    # Spin-2 fields
    f2 = [nmt.NmtField(n, [e1, e2], n_iter=0) for n, e1, e2 in zip(nmap, e1map, e2map)]

    # DD power spectra
    w_dd = {}
    for p in range(nbins):
        w_dd[p] = nmt.NmtWorkspace()
        w_dd[p].compute_coupling_matrix(f0[p],f0[p], b)
    cl_dd = np.array([[w_dd[p2].decouple_cell(nmt.compute_coupled_cell(f0[p1], f0[p2]))
                        for p2 in range(nbins)]
                    for p1 in range(nbins)])

    # w_dd = nmt.NmtWorkspace()
    # w_dd.compute_coupling_matrix(f0, f0, b)
    # cl_dd = np.array([w_dd.decouple_cell(nmt.compute_coupled_cell(f0[p1], f0[p2]))[0]
    #                   for p1, p2 in pairs])
        
    # DM power spectra
    w_dl = {}
    for p in range(nbins):
        w_dl[p] = nmt.NmtWorkspace()
        w_dl[p].compute_coupling_matrix(f0[p], f2[p], b)
    cl_dm = np.array([[w_dl[p2].decouple_cell(nmt.compute_coupled_cell(f0[p1], f2[p2]))
                       for p2 in range(nbins)]
                      for p1 in range(nbins)])

    # MM power spectra
    w_ll = {}
    for p1, p2 in pairs:
        w_ll[f'{p1}{p2}'] = nmt.NmtWorkspace()
        w_ll[f'{p1}{p2}'].compute_coupling_matrix(f2[p1], f2[p2], b)
    cl_mm = np.array([w_ll[f'{p1}{p2}'].decouple_cell(nmt.compute_coupled_cell(f2[p1], f2[p2]))
                       for p1, p2 in pairs])
        
    #d_values = np.array([hp.anafast(np.asarray([dmap[p1],e1map[p1],e2map[p1]]),np.asarray([dmap[p2],e1map[p2],e2map[p2]]), pol=True) for p1,p2 in pairs])
    #
    #cl_md_d = np.copy(d_values[:,3])
    #
    #for i, (p1,p2) in enumerate(pairs):
    #    if p1 != p2:
    #        cl_md_d[i] = np.array(hp.anafast(np.asarray([dmap[p2],e1map[p2],e2map[p2]]), np.asarray([dmap[p1],e1map[p1],e2map[p1]])))[3]

    return shotnoise, pairs, nz_tot, z_nz, cl_dd, cl_dm, cl_mm, cl_dd_t, cl_dm_t, cl_md_t, cl_mm_t


def compute_data(sim_path,source=1, output_path=None, nside=128, max_files=None, downsampling=1, zbins=[-1,0.15,1], nz_h = 50, nz_min=0, nz_max=None, code=None):
    ''' Method to compute the values needed for CCL test plots.
    
    Args:
        sim_path (str): Path where the CoLoRe simulation is located.
        source (int, optional): Source of which to compute data (default: 1)
        output_path (str, optional): Output where to save the data (default: { sim_path }/ccl_data/{ datetime.now() }/)
    '''
    if code != 'namaster':
        raise ValueError('Configuration set to run namaster code.')

    id_ = datetime.today().strftime('%Y%m%d_%H%M%S')
    if not output_path:
        output_path = sim_path + f"/ccl_data/{ id_ }"
    
    os.makedirs(output_path, exist_ok = True)

    log.debug(f'Computing data for:\nsim_path: { sim_path }\nsource: { source }\noutput_path: { output_path }')

    shotnoise, pairs, nz_tot, z_nz, cl_dd, cl_dm, cl_mm, cl_dd_t, cl_dm_t, cl_md_t, cl_mm_t = compute_all_cls(sim_path, source, nside, max_files, downsampling, zbins, nz_h, nz_min, nz_max)

    cl_dd_d = np.array([cl_dd[0, 0, 0, :],
                        cl_dd[0, 1, 0, :],
                        cl_dd[1, 1, 0, :]])
    cl_mm_d = cl_mm[:, 0, :]
    cl_bb_d = cl_mm[:, 3, :]
    cl_dm_d = np.array([cl_dm[0, 0, 0, :],
                        cl_dm[0, 1, 0, :],
                        cl_dm[1, 1, 0, :]])
    cl_md_d = np.array([cl_dm[0, 0, 0, :],
                        cl_dm[1, 0, 0, :],
                        cl_dm[1, 1, 0, :]])
    cl_mb_d = cl_mm[:, 1, :]
    cl_db_d = np.array([cl_dm[0, 0, 1, :],
                        cl_dm[0, 1, 1, :],
                        cl_dm[1, 1, 1, :]])

    savetofile(output_path, (pairs, shotnoise, nz_tot, z_nz, cl_dd_d, cl_dd_t, cl_dm_d, cl_dm_t, cl_md_d, cl_md_t, cl_mm_d, cl_mm_t, cl_bb_d, cl_mb_d, cl_db_d), ('pairs', 'shotnoise', 'nz_tot', 'z_nz', 'cl_dd_d', 'cl_dd_t', 'cl_dm_d', 'cl_dm_t', 'cl_md_d', 'cl_md_t', 'cl_mm_d', 'cl_mm_t','cl_bb_d','cl_mb_d','cl_db_d'))

    info = {
        'id'            : id_,
        'source'        : source,
        'nside'         : nside,
        'max_files'     : max_files,
        'downsampling'  : downsampling,
        'zbins'         : zbins,
        'nz_h'          : nz_h,
        'nz_min'        : nz_min,
        'nz_max'        : nz_max,
        'code'          : 'namaster'
    }

    with open(output_path + '/INFO.json','w') as outfile:
        json.dump(info, outfile)


def main(): #pragma: no cover
    parser = argparse.ArgumentParser(description="Save values to compute CCL test into .dat files")
    parser.add_argument("-p","--path", required=True, type=str, help="Path of ColoRe run")
    parser.add_argument("-o","--output", required=False, type=str, default=None, help="Path for output files")
    
    parser.add_argument("--source",       required=False, type=int, default=1, help="Sources to be computed")
    parser.add_argument("--nside",        required=False, type=int , default=128 , help="nside to use ")
    parser.add_argument("--max_files",    required=False, type=int, default=None , help="number of srcs files to consider (default: None, consider all the files)")
    parser.add_argument("--downsampling", required=False, type=float , default=1 , help="downsampling to apply to the data")
    parser.add_argument("--zbins",        required=False, type=float, nargs='+', default=[0,0.15,1] , help="defines the binning in redshift of the analysis")
    parser.add_argument("--nz_h",         required=False, type=int , default=50 , help="pixelization of the redshift analysis ")
    parser.add_argument("--nz_min",       required=False, type=float , default=0 , help="min redshift for the redshfit analysis")
    parser.add_argument("--nz_max",       required=False, type=float , default=None , help="max redshift for the redshfit analysis")

    args = parser.parse_args()
    
    path   = args.path
    output  = args.output

    compute_data(path, args.source, output, args.nside, args.max_files, args.downsampling, args.zbins, args.nz_h, args.nz_min, args.nz_max )

if __name__ == '__main__': #pragma: no cover
    main()