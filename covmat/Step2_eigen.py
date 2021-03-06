from numpy import *
import glob
from astropy.io import fits
import subprocess

def save_img(dat, imname):
    subprocess.getoutput("rm -f " + imname)
    fitsobj = fits.HDUList()
    hdu = fits.PrimaryHDU()
    hdu.data = dat
    fitsobj.append(hdu)
    fitsobj.writeto(imname)
    fitsobj.close()


for include_pecvel in [0, 1]:
    the_count = 0

    d_mBx1c_dsys = []
    nsne = 740

    for fl in glob.glob("C*fits"):
        if (fl.count("stat") == 0) or fl == "C_stat_new.fits":
            if (fl.count("pecvel") == 0) or include_pecvel:
                print("include_pecvel", include_pecvel, fl)
                f = fits.open(fl)
                dat = f[0].data
                f.close()

                evals, evecs = linalg.eig(dat)
                evecs = transpose(evecs)

                evals = real(evals)
                evecs = real(evecs)

                inds = argsort(evals)[::-1]

                for ind in inds:
                    if evals[ind] > evals[inds[0]]*1e-6:
                        d_mBx1c_dsys.append([])
                        for i in range(nsne):
                            d_mBx1c_dsys[-1].append((evecs[ind]*sqrt(evals[ind]))[i*3:(i+1)*3])

    d_mBx1c_dsys = array(d_mBx1c_dsys)
    print(d_mBx1c_dsys.shape)

    save_img(d_mBx1c_dsys, "d_mBx1c_dsys_pecvel=%i.fits" % include_pecvel)


