from numpy import *
import pyfits
from DavidsNM import save_img
from FileRead import readcol

[NA, zcmb, zhel, dz, mb, dmb, x1, dx1, color, dcolor, mass, dmass, cov_m_s, cov_m_c, cov_s_c] = readcol("jla_lcparams.txt", 'a,fff,ff,ff,ff,ff,fff')

f = pyfits.open("orig/C_stat.fits")
dat = f[0].data
f.close()

[sigma_coh, sigma_lens, sigma_z] = readcol("sigma_mu.txt", 'fff')


sigma_pecvel = (5. * 150. / 3e5) / (log(10.) * sigma_z)

for i in range(len(zcmb)):
    dat[i*3,i*3] -= dmb[i]**2 - sigma_coh[i]**2. - sigma_lens[i]**2. - sigma_pecvel[i]**2.
    dat[i*3+1,i*3+1] -= dx1[i]**2
    dat[i*3+2,i*3+2] -= dcolor[i]**2

    dat[i*3,i*3+1] -= cov_m_s[i]
    dat[i*3+1,i*3] -= cov_m_s[i]

    dat[i*3,i*3+2] -= cov_m_c[i]
    dat[i*3+2,i*3] -= cov_m_c[i]

    dat[i*3+1,i*3+2] -= cov_s_c[i]
    dat[i*3+2,i*3+1] -= cov_s_c[i]


save_img(dat, "C_stat.fits")
