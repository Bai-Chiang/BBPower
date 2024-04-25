import numpy as np

#Foreground model
A_sync_BB = 2.0
EB_sync = 2.
alpha_sync_EE = -0.6
alpha_sync_BB = -0.4
beta_sync = -3.
nu0_sync = 23.

A_dust_BB = 5.0
EB_dust = 2.
alpha_dust_EE = -0.42
alpha_dust_BB = -0.2
beta_dust = 1.59
temp_dust = 19.6
nu0_dust = 353.

Alens = 1.0

band_names = ['LF1', 'LF2', 'MF1', 'MF2', 'UHF1', 'UHF2']


#CMB spectrum
def fcmb(nu):
    x = 0.017608676067552197*nu
    ex = np.exp(x)
    return ex*(x/(ex-1))**2


#All spectra
def comp_sed(nu,nu0,beta,temp,typ):
    if typ == 'cmb':
        return fcmb(nu)
    elif typ == 'dust':
        x_to=0.04799244662211351*nu/temp
        x_from=0.04799244662211351*nu0/temp
        return (nu/nu0)**(1+beta)*(np.exp(x_from)-1)/(np.exp(x_to)-1)*fcmb(nu0)
    elif typ == 'sync':
        return (nu/nu0)**beta*fcmb(nu0)
    return None


#Component power spectra
def dl_plaw(A,alpha,ls):
    return A*((ls+0.001)/80.)**alpha


def get_nmt_binning(nside, delta_ell):
    """
    Produce an NmtBin object given nside and the multipole separation.
    """
    import pymaster as nmt

    bin_low = np.arange(0, 3*nside, delta_ell)
    bin_high = bin_low + delta_ell - 1
    bin_high[-1] = 3*nside - 1

    return nmt.NmtBin.from_edges(bin_low, bin_high + 1)


def read_camb(fname, lmax):
    larr_all = np.arange(lmax+1)
    l,dtt,dee,dbb,dte = np.loadtxt(fname,unpack=True)
    l = l.astype(int)
    msk = l <= lmax
    l = l[msk]
    dltt = np.zeros(len(larr_all))
    dltt[l] = dtt[msk]
    dlee = np.zeros(len(larr_all))
    dlee[l] = dee[msk]
    dlbb = np.zeros(len(larr_all))
    dlbb[l] = dbb[msk]
    dlte = np.zeros(len(larr_all))
    dlte[l] = dte[msk]
    return dltt,dlee,dlbb,dlte


#Bandpasses 
class Bpass(object):
    def __init__(self, name, fname):
        self.name = name
        self.nu, self.bnu = np.loadtxt(fname, unpack=True)
        self.dnu = np.zeros_like(self.nu)
        self.dnu[1:] = np.diff(self.nu)
        self.dnu[0] = self.dnu[1]
        # CMB units
        norm = np.sum(self.dnu*self.bnu*self.nu**2*fcmb(self.nu))
        self.bnu /= norm

    def convolve_sed(self, f):
        sed = np.sum(self.dnu*self.bnu*self.nu**2*f(self.nu))
        return sed


def read_beam_window(fname, lmax):
    bls = np.loadtxt(fname)[:, 1]
    ls = np.loadtxt(fname)[:, 0]
    return bls[ls <= lmax]


def get_component_spectra(
        lmax, params,
        fname_camb_lens_nobb="./examples/data/camb_lens_nobb.dat"
    ):
    """
    """
    larr_all = np.arange(lmax + 1)
    dls_sync_ee = dl_plaw(params["synch"]["A_s_EE"],
                          params["synch"]["alpha_s_EE"],
                          larr_all)
    dls_sync_bb = dl_plaw(params["synch"]["A_s_BB"],
                          params["synch"]["alpha_s_BB"],
                          larr_all)
    dls_dust_ee = dl_plaw(params["dust"]["A_d_EE"],
                          params["dust"]["alpha_d_EE"],
                          larr_all)
    dls_dust_bb = dl_plaw(params["dust"]["A_d_BB"],
                          params["dust"]["alpha_d_BB"],
                          larr_all)
    _, dls_cmb_ee, dls_cmb_bb, _ = read_camb(fname_camb_lens_nobb, lmax)

    return (dls_sync_ee, dls_sync_bb,
            dls_dust_ee, dls_dust_bb,
            dls_cmb_ee, params["CMB"]["A_lens"]*dls_cmb_bb)


def get_convolved_seds(band_names, bpss, params):
    """
    """
    nfreqs = len(band_names)
    seds = np.zeros([3, nfreqs])
    for ib, n in enumerate(band_names):
        b = bpss[n]
        seds[0, ib] = b.convolve_sed(
            lambda nu : comp_sed(nu, None, None, None, 'cmb')
        )
        seds[1, ib] = b.convolve_sed(
            lambda nu : comp_sed(nu,
                                 params["synch"]["nu0_s"],
                                 params["synch"]["beta_s"],
                                 None, 'sync')
        )
        seds[2, ib] = b.convolve_sed(
            lambda nu : comp_sed(nu,
                                 params["dust"]["nu0_d"],
                                 params["dust"]["beta_d"],
                                 params["dust"]["T_d"],
                                 'dust')
        )
    return seds
