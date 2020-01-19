import pystan
import pickle as pickle
from numpy import *
from matplotlib import use
use("PDF")
import matplotlib.pyplot as plt
from astropy.io import fits
import argparse
from scipy.interpolate import interp1d
import time
import astropy.io.ascii as ascii
import subprocess


def radectoxyz(RAdeg, Decdeg):
    x = cos(Decdeg/(180./pi))*cos(RAdeg/(180./pi))
    y = cos(Decdeg/(180./pi))*sin(RAdeg/(180./pi))
    z = sin(Decdeg/(180./pi))

    return array([x, y, z], dtype=float64)



def get_dz(RAdeg, Decdeg):
    
    dzCMB = 371.e3/299792458.
    CMBcoordsRA = 168.01190437
    CMBcoordsDEC = -6.98296811
          

    CMBxyz = radectoxyz(CMBcoordsRA, CMBcoordsDEC)
    inputxyz = radectoxyz(RAdeg, Decdeg)
    
    dz = dzCMB*dot(CMBxyz, inputxyz)
    dv = dzCMB*dot(CMBxyz, inputxyz)*299792.458

    print("Add this to z_helio to lowest order:")
    print(dz, dv)

    return dz


def get_zCMB(RAdeg, Decdeg, z_helio):
    dz = -get_dz(RAdeg, Decdeg)

    one_plus_z_pec = sqrt((1. + dz)/(1. - dz))
    one_plus_z_CMB = (1 + z_helio)/one_plus_z_pec
    return one_plus_z_CMB - 1.

def get_dot_CMB(RAdeg, Decdeg):
    CMBcoordsRA = 168.01190437
    CMBcoordsDEC = -6.98296811

    CMBxyz = radectoxyz(CMBcoordsRA, CMBcoordsDEC)
    inputxyz = radectoxyz(RAdeg, Decdeg)

    return dot(CMBxyz, inputxyz)


def save_img(dat, imname):
    subprocess.getoutput("rm -f " + imname)
    fitsobj = fits.HDUList()
    hdu = fits.PrimaryHDU()
    hdu.data = dat
    fitsobj.append(hdu)
    fitsobj.writeto(imname)
    fitsobj.close()


def get_redshifts(redshifts):
    appended_redshifts = arange(0., 2.51, 0.1)
    tmp_redshifts = concatenate((redshifts, appended_redshifts))
    
    sort_inds = list(argsort(tmp_redshifts))
    unsort_inds = [sort_inds.index(i) for i in range(len(tmp_redshifts))]
    
    tmp_redshifts = sort(tmp_redshifts)
    redshifts_sort_fill = sort(concatenate((tmp_redshifts, 0.5*(tmp_redshifts[1:] + tmp_redshifts[:-1]))))
    
    return redshifts, redshifts_sort_fill, unsort_inds, len(appended_redshifts)


def get_redshift_coeffs(zcmb, SNset, popmodel):
    
    if popmodel == 0:
        redshift_coeffs = ones([len(zcmb), 3], dtype=float64)
        redshift_coeffs[:, 1] = zcmb
        redshift_coeffs[:, 2] = zcmb - zcmb**2. # Slightly decorrelate
        
    if popmodel == 1:
        return ones([len(zcmb), 1], dtype=float64)

    if popmodel == 2:
        redshift_coeffs = zeros([len(zcmb), len(unique(SNset))], dtype=float64)
        for i, id in enumerate(unique(SNset)):
            redshift_coeffs[:,i] = (SNset == id)

    if popmodel > 2:
        npersample = popmodel - 1

        redshift_coeffs = zeros([len(zcmb), len(unique(SNset))*npersample - 1], dtype=float64)
        assert sum(SNset == 4) < 30 # Just checking that SNset 4 is HST SNe

        the_pos = 0
        for id in unique(SNset):
            minz = (zcmb[where(SNset == id)]).min()
            maxz = (zcmb[where(SNset == id)]).max()
            dz = maxz - minz

            if id < 4:
                for j in range(npersample):
                    yvals = zeros(npersample, dtype=float64)
                    yvals[j] = 1.

                    ifn = interp1d(linspace(minz - 1e-8, maxz + 1e-8, npersample), yvals, kind = 'linear', fill_value = 0, bounds_error = False)
                    redshift_coeffs[:,the_pos] = (SNset == id)*ifn(zcmb)
                    the_pos += 1
            else:
                redshift_coeffs[:,the_pos] = (SNset == id)
                the_pos += 1

    return redshift_coeffs

        
def initfn():
    if args.cosmomodel == 1:
        Ominit = 0.3 + random.random()*0.1
        OLinit = 0.7 + random.random()*0.1
        q0init = 0.
        j0init = 0.
    if args.cosmomodel == 2:
        Ominit = 0.3 + random.random()*0.1
        OLinit = 0.01
        q0init = 0.01
        j0init = 0.01
    if args.cosmomodel == 3 or args.cosmomodel == 5:
        Ominit = 0.01
        OLinit = 0.01
        q0init = -0.5 + random.random()*0.1
        j0init = random.normal()*0.1
    if args.cosmomodel == 4:
        Ominit = 0.3 + random.random()*0.1
        OLinit = -1. + random.random()*0.1
        q0init = 0.01
        j0init = 0.01
    
        

    return dict(alpha = 0.12 + random.normal()*0.01,
                beta = 3. + random.normal()*0.1,
                delta = random.normal()*0.01,
                
                Om = Ominit,
                OL = OLinit,
                q0 = q0init,
                j0 = j0init,
                q0m = -0.5 + random.random()*0.1,
                q0d = random.normal(),

                sigma_M0 = random.random()*0.01 + 0.1,
                sigma_x10 = random.random()*0.1 + 1.,
                sigma_c0 = random.random()*0.01 + 0.06,
                
                M0 = -19.1 + random.random()*0.1,
                x10 = random.normal(size = len(redshift_coeffs[0]))*0.1,
                c0 = random.normal(size = len(redshift_coeffs[0]))*0.01,
                
                true_x1 = random.normal(size = nsne)*0.1,
                true_c = random.normal(size = nsne)*0.01,
                
                calibs = random.normal(size = d_mBx1c_dsys.shape[2])*0.1)
                


parser = argparse.ArgumentParser()
parser.add_argument("--cosmomodel", type=int, help="1 = Om/OL, 2 = FlatLCDM, 3 = q0/j0, 4 = q0m/q0d/j0")
parser.add_argument("--popmodel", type=int, help="0 = z, z^2, 1 = const, 2 = const by sample, 3 = linear by sample")
parser.add_argument("--hostmass", type=int, help="host mass? 1 = yes, 0 = no")
parser.add_argument("--includepecvelcov", type=int, help="include peculiar velocity covariance matrix? 1 = yes, 0 = no")
parser.add_argument("--ztype", type=str, help="redshift type to use for comoving distance. zcmbpecvel, zcmb, or zhelio")
parser.add_argument("--nMCMCchains", type=int, help="number of chains to run")
parser.add_argument("--nMCMCsamples", type=int, help="number of samples per chain; first half is discarded")
parser.add_argument("--min_Om", type=float, help="minimum Omega_m", default = 0)
parser.add_argument("--saveperSN", type=int, help="Save per-SN parameters in pickle file?", default = 1)
parser.add_argument("--savestan", type=int, help="Save Stan data in pickle", default = 1)


args = parser.parse_args()
print("args ", args)


lcparams = ascii.read("../covmat/jla_lcparams.txt")
sigmamu = ascii.read("../covmat/sigma_mu.txt", names = ["sigma_coh", "sigma_lens", "z"])


assert all(abs(sigmamu["z"] - lcparams["zcmb"]) < 0.02)

dmb = sqrt(lcparams["dmb"]**2. - sigmamu["sigma_coh"]**2.)

plt.plot(lcparams["zcmb"], dmb, '.')
plt.savefig("dmb_vs_z.pdf")
plt.close()

f = fits.open("../covmat/d_mBx1c_dsys_pecvel=%i.fits" % args.includepecvelcov)
d_mBx1c_dsys = f[0].data
f.close()

d_mBx1c_dsys = transpose(d_mBx1c_dsys, axes = [1, 2, 0])


dot_CMB = array([get_dot_CMB(lcparams["ra"][i], lcparams["dec"][i]) for i in range(len(lcparams["ra"]))])

all_z = dict(zcmbpecvel = lcparams["zcmb"],
             zcmb = array([get_zCMB(lcparams["ra"][i], lcparams["dec"][i], lcparams["zhel"][i]) for i in range(len(lcparams["ra"]))]),
             zhelio = lcparams["zhel"])
assert args.ztype in all_z, "available z keys: " + str(all_z.keys())


print(d_mBx1c_dsys.shape)


nsne = len(lcparams["zcmb"])

obs_mBx1c = zeros([nsne, 3], dtype=float64)
obs_mBx1c_cov = zeros([nsne, 3,3], dtype=float64)

for i in range(nsne):
    obs_mBx1c[i] = [lcparams["mb"][i], lcparams["x1"][i], lcparams["color"][i]]
    obs_mBx1c_cov[i] = [[dmb[i]**2., lcparams["cov_m_s"][i], lcparams["cov_m_c"][i]],
                        [lcparams["cov_m_s"][i], lcparams["dx1"][i]**2., lcparams["cov_s_c"][i]],
                        [lcparams["cov_m_c"][i], lcparams["cov_s_c"][i], lcparams["dcolor"][i]**2.]]

save_img(obs_mBx1c_cov, "obs_mBx1c_cov.fits")

redshifts, redshifts_sort_fill, unsort_inds, nzadd = get_redshifts(all_z[args.ztype]) # CMB for this one, helio for the other one!

redshift_coeffs = get_redshift_coeffs(all_z[args.ztype], lcparams["set"], args.popmodel)

for i in range(len(redshift_coeffs[0])):
    plt.plot(lcparams["zcmb"], redshift_coeffs[:,i] + random.normal(size = nsne)*0.01, '.', label = str(i))

plt.ylim(-0.2, 1.2)
plt.legend(loc = 'best')
plt.xscale('log')
plt.savefig("redshift_coeffs_%i.pdf" % args.popmodel)
plt.close()



stan_data = dict(n_sne = nsne, n_calib = d_mBx1c_dsys.shape[2], nzadd = nzadd, n_x1c_star = len(redshift_coeffs[0]),
                 zhelio = lcparams["zhel"], zcmb = all_z[args.ztype], dot_CMB = dot_CMB, redshifts_sort_fill = redshifts_sort_fill, unsort_inds = unsort_inds,
                 redshift_coeffs = redshift_coeffs,
                 obs_mBx1c = obs_mBx1c, obs_mBx1c_cov = obs_mBx1c_cov,
                 d_mBx1c_d_calib = d_mBx1c_dsys,
                 obs_mass = lcparams["3rdvar"], obs_dmass = lcparams["d3rdvar"],
                 cosmomodel = args.cosmomodel, min_Om = args.min_Om, host_mass_relation = args.hostmass)

plt.subplot(2,1,1)
plt.hist(lcparams["3rdvar"])
plt.subplot(2,1,2)
plt.hist(lcparams["d3rdvar"], bins = 20)
plt.savefig("mass.pdf")
plt.close()

print("Ready to sample", time.asctime())

fit = pystan.stan(file = "../stan_code.txt", data=stan_data,
                  iter=args.nMCMCsamples, chains=args.nMCMCchains, n_jobs = args.nMCMCchains, refresh = int(min(100, args.nMCMCsamples/20)), init = initfn)

print("Done with sampling", time.asctime())
print(fit)
print("Done with printing", time.asctime())

fit_params = fit.extract(permuted = True)

print("Done with extracting", time.asctime())

if args.saveperSN:
    pass
else:
    for key in fit_params:
        if fit_params[key].size > 100 * args.nMCMCsamples * args.nMCMCchains:
            print("Deleting ", key)
            fit_params[key] = array([], dtype=float64)


if args.savestan:
    pass
else:
    stan_data = {}


pickle.dump((stan_data, fit_params), open("results.pickle", 'wb'))
print("Done!", time.asctime())
