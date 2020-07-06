import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
from astropy.io import fits
import healpy as hp
import sys
import os


def read_ascii(fname):
    ifile = 0

    d = np.loadtxt(fname + "_%d.txt" % ifile, unpack=True)
    has_shear = False
    if len(d) == 7:
        has_shear = True
        type_arr, ra_arr, dec_arr, z0_arr, rsd_arr, e1_arr, e2_arr = d
    else:
        type_arr, ra_arr, dec_arr, z0_arr, rsd_arr = d
        e1_arr = None
        e2_arr = None

    ifile = 1
    while os.path.isfile(fname + "_%d.h5" % ifile):
        d = np.loadtxt(fname + "_%d.txt" % ifile, unpack=True)
        if has_shear:
            typ, ra, dec, z0, rsd, e1, e2 = d
        else:
            typ, ra, dec, z0, rsd = d
        type_arr = np.concatenate((type_arr, type_arr))
        ra_arr = np.concatenate((ra_arr, ra_arr))
        dec_arr = np.concatenate((dec_arr, dec_arr))
        z0_arr = np.concatenate((z0_arr, z0_arr))
        rsd_arr = np.concatenate((rsd_arr, rsd_arr))
        if has_shear:
            e1_arr = np.concatenate((e1_arr, e1_arr))
            e2_arr = np.concatenate((e2_arr, e2_arr))
        ifile += 1

    return ra_arr, dec_arr, z0_arr, rsd_arr, \
        type_arr, e1_arr, e2_arr, None, None, None


def read_hdf5(fname, ipop):
    ifile = 0

    ff = h5.File(fname + "_%d.h5" % ifile, "r")
    tabname = '/sources%d' % ipop

    ra_arr = ff[tabname]['RA']
    dec_arr = ff[tabname]['DEC']
    z0_arr = ff[tabname]['Z_COSMO']
    rsd_arr = ff[tabname]['DZ_RSD']

    has_shear = 'E1' in ff[tabname].keys()
    if has_shear:
        e1_arr = ff[tabname]['E1']
        e2_arr = ff[tabname]['E2']
    else:
        e1_arr = None
        e2_arr = None

    ifile = 1
    while os.path.isfile(fname + "_%d.h5" % ifile):
        ff = h5.File(fname + "_%d.h5" % ifile, "r")
        tabname = '/sources%d' % ipop

        ra_arr = np.concatenate((ra_arr, ff[tabname]['RA']))
        dec_arr = np.concatenate((dec_arr, ff[tabname]['DEC']))
        z0_arr = np.concatenate((z0_arr, ff[tabname]['Z_COSMO']))
        rsd_arr = np.concatenate((rsd_arr, ff[tabname]['DZ_RSD']))
        if has_shear:
            e1_arr = np.concatenate((e1_arr, ff[tabname]['E1']))
            e2_arr = np.concatenate((e2_arr, ff[tabname]['E2']))
        ifile += 1

    type_arr = ipop*np.ones(len(ra_arr))

    return ra_arr, dec_arr, z0_arr, rsd_arr, \
        type_arr, e1_arr, e2_arr, None, None, None


def read_fits(fname):
    ifile = 0

    fi = fits.open(fname + "_%d.fits" % ifile)
    data = (fi[1]).data
    with_skewers = len(fi) > 2

    if with_skewers:
        dens_skewers = (fi[2]).data
        vrad_skewers = (fi[3]).data
        data_skewers = (fi[4]).data
    else:
        dens_skewers = None
        vrad_skewers = None
        data_skewers = None

    ra_arr = data['RA']
    dec_arr = data['DEC']
    z0_arr = data['Z_COSMO']
    rsd_arr = data['DZ_RSD']
    has_shear = 'E1' in data.names
    if has_shear:
        e1_arr = data['E1']
        e2_arr = data['E2']
    else:
        e1_arr = None
        e2_arr = None
    type_arr = data['TYPE']

    ifile = 1
    while os.path.isfile(fname + "_%d.fits" % ifile):
        fi = fits.open(fname + "_%d.fits" % ifile)
        data = (fi[1]).data
        if with_skewers:
            dens_skw = (fi[2]).data
            vrad_skw = (fi[3]).data

        ra_arr = np.concatenate((ra_arr, data['RA']))
        dec_arr = np.concatenate((dec_arr, data['DEC']))
        z0_arr = np.concatenate((z0_arr, data['Z_COSMO']))
        rsd_arr = np.concatenate((rsd_arr, data['DZ_RSD']))
        if has_shear:
            e1_arr = np.concatenate((e1_arr, data['E1']))
            e2_arr = np.concatenate((e2_arr, data['E2']))
        type_arr = np.concatenate((type_arr, data['TYPE']))
        if with_skewers:
            dens_skewers = np.concatenate((dens_skewers, dens_skw))
            vrad_skewers = np.concatenate((vrad_skewers, vrad_skw))
        ifile += 1

    return ra_arr, dec_arr, z0_arr, rsd_arr, type_arr, \
        e1_arr, e2_arr, dens_skewers, vrad_skewers, data_skewers


if len(sys.argv) != 3:
    print("Usage: read_colore.py file_name file_format "
          "('ASCII', 'FITS' or 'HDF5')")
    exit(1)
fname = sys.argv[1]
fmt = sys.argv[2]

if fmt == 'ASCII':
    ra_arr, dec_arr, z_arr, rsd_arr, type_arr, \
        e1_arr, e2_arr, dskw_arr, vskw_arr, data_skw = read_ascii(fname)
elif fmt == 'FITS':
    ra_arr, dec_arr, z_arr, rsd_arr, type_arr, \
        e1_arr, e2_arr, dskw_arr, vskw_arr, data_skw = read_fits(fname)
elif fmt == 'HDF5':
    ra_arr, dec_arr, z_arr, rsd_arr, type_arr, \
        e1_arr, e2_arr, dskw_arr, vskw_arr, data_skw = read_hdf5(fname, 1)

nside = 64
npix = hp.nside2npix(nside)
ipixs = hp.ang2pix(nside, np.radians(90-dec_arr), np.radians(ra_arr))
mp = np.histogram(ipixs, bins=npix, range=[0, npix])[0]
hp.mollview(mp, title='$N_g(\\hat{\\bf n})$')
mpr = np.histogram(ipixs, bins=npix, range=[0, npix], weights=rsd_arr)[0]
hp.mollview(mpr/mp, title='$v_r$')
if e1_arr is not None:
    mpe1 = np.histogram(ipixs, bins=npix, range=[0, npix], weights=e1_arr)[0]
    mpe2 = np.histogram(ipixs, bins=npix, range=[0, npix], weights=e2_arr)[0]
    hp.mollview(mpe1/mp, title='$e_1$')
    hp.mollview(mpe2/mp, title='$e_2$')
plt.figure()
plt.hist(z_arr, bins=100, histtype='step')
plt.xlabel('$z$', fontsize=16)
plt.ylabel('$N(z)$', fontsize=16)

if dskw_arr is not None:
    # Plot skewers
    nside_skw = 64
    ipix_skw = 1000
    id_in_pix = np.where(hp.ang2pix(nside_skw,
                                    np.radians(90 - dec_arr),
                                    np.radians(ra_arr))
                         == ipix_skw)[0]
    cols = ['r', 'g', 'b', 'y', 'k']
    plt.figure()
    plt.title('Skewers')
    plt.xlabel('$z$')
    plt.ylabel('$\\delta$')
    for i in id_in_pix:
        plt.plot(data_skw['Z'], dskw_arr[i, :], cols[i % 5]+'-')
        plt.plot([z_arr[i], z_arr[i]],
                 [1.1 * np.amin(dskw_arr[i, :]),
                  1.1 * np.amax(dskw_arr[i, :])],
                 cols[i % 5]+'--')
    plt.figure()
    plt.title('Skewers')
    plt.xlabel('$z$')
    plt.ylabel('$v_r$')
    for i in id_in_pix:
        plt.plot(data_skw['Z'], vskw_arr[i, :], cols[i % 5]+'-')
        plt.plot([z_arr[i], z_arr[i]],
                 [1.1 * np.amin(vskw_arr[i, :]),
                  1.1 * np.amax(vskw_arr[i, :])],
                 cols[i % 5]+'--')
plt.show()
