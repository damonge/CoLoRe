import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import pyccl as ccl
from scipy.integrate import simps


h = 0.7
cosmo = ccl.Cosmology(Omega_c=0.25, Omega_b=0.05, h=h, n_s=0.96, sigma8=0.8)

for i in range(1,3):
    z_w, w = np.loadtxt(f"kz{i}.txt", unpack=True)
    chi_w = ccl.comoving_radial_distance(cosmo, 1./(1+z_w))
    h_w = ccl.h_over_h0(cosmo, 1./(1+z_w))*h/2997.92458
    mp = hp.read_map(f"out_custom_s{i}.fits", verbose=False)
    t = ccl.Tracer()
    t.add_tracer(cosmo, kernel=(chi_w, w*h_w))

    zpk = np.array([0.0, 0.1, 0.2, 0.3, 0.4])
    pks = []
    for z in zpk:
        fname = "out_pk_custom_pop%d_z%.3lf.txt" % (i-1, z)
        k, pk, _, pmm = np.loadtxt(fname, unpack=True)
        pk[pmm<1E-30] = 0
        pk[pk<0] = 0
        pk = pk[(k>=3E-4) & (k<=3E-1)]
        k = k[(k>=3E-4) & (k<=3E-1)]
        pks.append(pk)
    pks = np.array(pks)[::-1]
    zpk = zpk[::-1]
    apk = 1./(1+zpk)
    pk = ccl.Pk2D(a_arr=apk, lk_arr=np.log(k*h),
                  pk_arr=np.log(pks/h**3), is_logp=True,
                  extrap_order_hik=1, extrap_order_lok=1)
    assert pk.eval(1, 1, cosmo) > 0

    mp = hp.read_map(f"out_custom_s{i}.fits", verbose=False)
    nside = hp.npix2nside(len(mp))
    ls = np.arange(3*nside)
    cld = hp.anafast(mp)
    clt = ccl.angular_cl(cosmo, t, t, ls, p_of_k_a=pk)

    plt.figure()
    plt.errorbar(ls, cld, yerr=cld/np.sqrt(ls+0.5), fmt='k.')
    plt.plot(ls, clt, 'r-')
    #plt.xlim([2, nside])
    plt.loglog()
plt.show()
