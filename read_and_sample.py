import pystan
import pickle as pickle
from numpy import *
from matplotlib import use
use("PDF")
import matplotlib.pyplot as plt
import pyfits
import sys
from scipy.interpolate import interp1d
import time
import astropy.io.ascii as ascii
import subprocess

def save_img(dat, imname):
    subprocess.getoutput("rm -f " + imname)
    fitsobj = pyfits.HDUList()
    hdu = pyfits.PrimaryHDU()
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
    if cosmomodel == 1:
        Ominit = 0.3 + random.random()*0.1
        OLinit = 0.7 + random.random()*0.1
        q0init = 0.
        j0init = 0.
    if cosmomodel == 2:
        Ominit = 0.3 + random.random()*0.1
        OLinit = 0.01
        q0init = 0.01
        j0init = 0.01
    if cosmomodel == 3:
        Ominit = 0.01
        OLinit = 0.01
        q0init = -0.5 + random.random()*0.1
        j0init = random.normal()*0.1
    if cosmomodel == 4:
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

                sigma_M0 = random.random()*0.01 + 0.1,
                sigma_x10 = random.random()*0.1 + 1.,
                sigma_c0 = random.random()*0.01 + 0.06,
                
                M0 = -19.1 + random.random()*0.1,
                x10 = random.normal(size = len(redshift_coeffs[0]))*0.1,
                c0 = random.normal(size = len(redshift_coeffs[0]))*0.01,
                
                true_x1 = random.normal(size = nsne)*0.1,
                true_c = random.normal(size = nsne)*0.01,
                
                calibs = random.normal(size = d_mBx1c_dsys.shape[2])*0.1)
                


cosmomodel = int(sys.argv[1]) # 1 = Om/OL, 2 = FlatLCDM
popmodel = int(sys.argv[2]) # 1 = const, 2 = const by sample, 3 = linear by sample
hostmass = int(sys.argv[3])
nMCMCchains = int(sys.argv[4])
nMCMCsamples = int(sys.argv[5])

if len(sys.argv) == 7:
    min_Om = float(sys.argv[6])
else:
    min_Om = 0.

lcparams = ascii.read("../covmat/jla_lcparams.txt")
sigmamu = ascii.read("../covmat/sigma_mu.txt", names = ["sigma_coh", "sigma_lens", "z"])


assert all(abs(sigmamu["z"] - lcparams["zcmb"]) < 0.02)

dmb = sqrt(lcparams["dmb"]**2. - sigmamu["sigma_coh"]**2.)

plt.plot(lcparams["zcmb"], dmb, '.')
plt.savefig("dmb_vs_z.pdf")
plt.close()

f = pyfits.open("../covmat/d_mBx1c_dsys.fits")
d_mBx1c_dsys = f[0].data
f.close()

d_mBx1c_dsys = transpose(d_mBx1c_dsys, axes = [1, 2, 0])



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

redshifts, redshifts_sort_fill, unsort_inds, nzadd = get_redshifts(lcparams["zcmb"]) # CMB for this one, helio for the other one!

redshift_coeffs = get_redshift_coeffs(lcparams["zcmb"], lcparams["set"], popmodel)

for i in range(len(redshift_coeffs[0])):
    plt.plot(lcparams["zcmb"], redshift_coeffs[:,i] + random.normal(size = nsne)*0.01, '.', label = str(i))

plt.ylim(-0.2, 1.2)
plt.legend(loc = 'best')
plt.xscale('log')
plt.savefig("redshift_coeffs_%i.pdf" % popmodel)
plt.close()



stan_data = dict(n_sne = nsne, n_calib = d_mBx1c_dsys.shape[2], nzadd = nzadd, n_x1c_star = len(redshift_coeffs[0]),
                 zhelio = lcparams["zhel"], zcmb = lcparams["zcmb"], redshifts_sort_fill = redshifts_sort_fill, unsort_inds = unsort_inds,
                 redshift_coeffs = redshift_coeffs,
                 obs_mBx1c = obs_mBx1c, obs_mBx1c_cov = obs_mBx1c_cov,
                 d_mBx1c_d_calib = d_mBx1c_dsys,
                 obs_mass = lcparams["3rdvar"], obs_dmass = lcparams["d3rdvar"],
                 cosmomodel = cosmomodel, min_Om = min_Om, host_mass_relation = hostmass)

plt.subplot(2,1,1)
plt.hist(lcparams["3rdvar"])
plt.subplot(2,1,2)
plt.hist(lcparams["d3rdvar"], bins = 20)
plt.savefig("mass.pdf")
plt.close()

print("Ready to sample", time.asctime())

fit = pystan.stan(file = "../stan_code.txt", data=stan_data,
                  iter=nMCMCsamples, chains=nMCMCchains, n_jobs = nMCMCchains, refresh = min(100, nMCMCsamples/20), init = initfn)

print("Done with sampling", time.asctime())
print(fit)
print("Done with printing", time.asctime())

fit_params = fit.extract(permuted = True)

print("Done with extracting", time.asctime())

pickle.dump((stan_data, fit_params), open("results.pickle", 'wb'))
print("Done!", time.asctime())
