from numpy import *
import glob
import pyfits
import commands

def save_img(dat, imname):
    commands.getoutput("rm -f " + imname)
    fitsobj = pyfits.HDUList()
    hdu = pyfits.PrimaryHDU()
    hdu.data = dat
    fitsobj.append(hdu)
    fitsobj.writeto(imname)
    fitsobj.close()

the_count = 0

d_mBx1c_dsys = []
nsne = 740

for fl in glob.glob("C*fits"):
    if (fl.count("stat") == 0) or fl == "C_stat_new.fits":
        print fl
        f = pyfits.open(fl)
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
print d_mBx1c_dsys.shape

save_img(d_mBx1c_dsys, "d_mBx1c_dsys.fits")


