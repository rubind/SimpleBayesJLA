from numpy import *
import pyfits
import astropy.io.ascii as ascii
import commands

def save_img(dat, imname):
    commands.getoutput("rm -f " + imname)
    fitsobj = pyfits.HDUList()
    hdu = pyfits.PrimaryHDU()
    hdu.data = dat
    fitsobj.append(hdu)
    fitsobj.writeto(imname)
    fitsobj.close()



lcparams = ascii.read("jla_lcparams.txt")
sigmamu = ascii.read("sigma_mu.txt", names = ["sigma_coh", "sigma_lens", "z"])



f = pyfits.open("C_stat.fits")
dat = f[0].data
f.close()



sigma_pecvel = (5. * 150. / 3e5) / (log(10.) * sigmamu["z"])

for i in range(len(lcparams["zcmb"])):
    dat[i*3,i*3] -= lcparams["dmb"][i]**2 - sigmamu["sigma_coh"][i]**2. - sigmamu["sigma_lens"][i]**2. - sigma_pecvel[i]**2.
    dat[i*3+1,i*3+1] -= lcparams["dx1"][i]**2
    dat[i*3+2,i*3+2] -= lcparams["dcolor"][i]**2

    dat[i*3,i*3+1] -= lcparams["cov_m_s"][i]
    dat[i*3+1,i*3] -= lcparams["cov_m_s"][i]

    dat[i*3,i*3+2] -= lcparams["cov_m_c"][i]
    dat[i*3+2,i*3] -= lcparams["cov_m_c"][i]

    dat[i*3+1,i*3+2] -= lcparams["cov_s_c"][i]
    dat[i*3+2,i*3+1] -= lcparams["cov_s_c"][i]


save_img(dat, "C_stat_new.fits")
