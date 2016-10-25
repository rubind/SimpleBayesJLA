import pystan
import cPickle as pickle
from numpy import *
import matplotlib.pyplot as plt
from FileRead import readcol
import pyfits
from DavidsNM import save_img
import sys

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

    if popmodel == 3:
        redshift_coeffs = zeros([len(zcmb), len(unique(SNset))*2 - 1], dtype=float64)
        assert sum(SNset == 4) < 30 # Just checkin that SNset 4 is HST SNe

        for i, id in enumerate(unique(SNset)):
            minz = (zcmb[where(SNset == id)]).min()
            maxz = (zcmb[where(SNset == id)]).max()
            dz = maxz - minz

            if id < 4:
                redshift_coeffs[:,2*i] = (SNset == id)*(maxz - zcmb)/dz
                redshift_coeffs[:,2*i+1] = (SNset == id)*(zcmb - minz)/dz
            else:
                redshift_coeffs[:,2*i] = (SNset == id)
    return redshift_coeffs

        


cosmomodel = int(sys.argv[1]) # 1 = Om/OL, 2 = FlatLCDM
popmodel = int(sys.argv[2]) # 1 = const, 2 = const by sample, 3 = linear by sample
nMCMCchains = int(sys.argv[3])
nMCMCsamples = int(sys.argv[4])


[NA, zcmb, zhel, dz, mb, dmb, x1, dx1, color, dcolor, mass, dmass, cov_m_s, cov_m_c, cov_s_c, SNset] = readcol("../covmat/jla_lcparams.txt", 'a,fff,ff,ff,ff,ff,fff,i')
[sigma_coh, sigma_lens, sigma_z] = readcol("../covmat/sigma_mu.txt", 'fff')

assert all(abs(sigma_z - zcmb) < 0.02)

dmb = sqrt(dmb**2. - sigma_coh**2.)

plt.plot(zcmb, dmb, '.')
plt.savefig("tmp.pdf")
plt.close()

f = pyfits.open("../covmat/d_mBx1c_dsys.fits")
d_mBx1c_dsys = f[0].data
f.close()

d_mBx1c_dsys = transpose(d_mBx1c_dsys, axes = [1, 2, 0])



print d_mBx1c_dsys.shape


nsne = len(zcmb)

obs_mBx1c = zeros([nsne, 3], dtype=float64)
obs_mBx1c_cov = zeros([nsne, 3,3], dtype=float64)

for i in range(nsne):
    obs_mBx1c[i] = [mb[i], x1[i], color[i]]
    obs_mBx1c_cov[i] = [[dmb[i]**2., cov_m_s[i], cov_m_c[i]],
                        [cov_m_s[i], dx1[i]**2., cov_s_c[i]],
                        [cov_m_c[i], cov_s_c[i], dcolor[i]**2.]]

save_img(obs_mBx1c_cov, "obs_mBx1c_cov.fits")

redshifts, redshifts_sort_fill, unsort_inds, nzadd = get_redshifts(zcmb) # CMB for this one, helio for the other one!

redshift_coeffs = get_redshift_coeffs(zcmb, SNset, popmodel)

for i in range(len(redshift_coeffs[0])):
    plt.plot(zcmb, redshift_coeffs[:,i] + random.normal(size = nsne)*0.01, '.', label = str(i))

plt.ylim(-0.2, 1.2)
plt.legend(loc = 'best')
plt.xscale('log')
plt.savefig("redshift_coeffs_%i.pdf" % popmodel)
plt.close()

stan_data = dict(n_sne = nsne, n_calib = d_mBx1c_dsys.shape[2], nzadd = nzadd, n_x1c_star = len(redshift_coeffs[0]),
                 zhelio = zhel, redshifts_sort_fill = redshifts_sort_fill, unsort_inds = unsort_inds,
                 redshift_coeffs = redshift_coeffs,
                 obs_mBx1c = obs_mBx1c, obs_mBx1c_cov = obs_mBx1c_cov,
                 d_mBx1c_d_calib = d_mBx1c_dsys,
                 obs_mass = mass, obs_dmass = dmass,
                 cosmomodel = cosmomodel)

plt.subplot(2,1,1)
plt.hist(mass)
plt.subplot(2,1,2)
plt.hist(dmass, bins = 20)
plt.savefig("mass.pdf")
plt.close()

fit = pystan.stan(file = "../stan_code.txt", data=stan_data,
                  iter=nMCMCsamples, chains=nMCMCchains, n_jobs = nMCMCchains, refresh = 100)


print fit

fit_params = fit.extract(permuted = True)

pickle.dump((stan_data, fit_params), open("results.pickle", 'wb'))
print "Done!"
