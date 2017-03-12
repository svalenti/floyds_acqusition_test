#!/usr/env python

from astropy.io import fits as pyfits
from astropy import wcs as pywcs
import numpy as np
from pyraf import iraf
import glob,os,re
import sys
import pylab as pl
import floacqastrodef
import scipy.signal as signal
from optparse import OptionParser

pl.ion()
fig = pl.figure(1,figsize=(6,6))
pl.rcParams['figure.figsize'] = 6, 6
pl.subplots_adjust(hspace = .01)
pl.subplots_adjust(wspace = .01)
pl.subplots_adjust(right=0.93, top=.93,left=.07, bottom=.07)
from matplotlib.image import imsave
pl.ion()
#   wcs command
#   pix1 = wcs1.wcs_sky2pix(sky0, 1)
#   sky0 = wcs.wcs_pix2sky(pixref, 1)

imglist = glob.glob('*fits')
for img in imglist:
    ################    update slit position 
    f = pyfits.open(img, mode='update')
    if not 'SPXX' in f[0].header:
        # making changes in data and/or header
        spxx = f[0].header['CRPIX1']
        spyy = f[0].header['CRPIX2']
        f[0].header['SPYY'] = spyy
        f[0].header['SPXX'] = spxx
        f.flush()  # changes are written back to original.fits
        f.close()  # closing the file will also flush an
    else:
        print 'slit position found'
        spxx = f[0].header['SPXX']
        spyy = f[0].header['SPYY']
        f.close()
    ################    
    # read it again 
    data, hdr0 = pyfits.getdata(img, header=True)

    if not 'xpos' in hdr0:
        f = pyfits.open(img, mode='update')
        # making changes in data and/or header

        # chose hd depending on the slit width
        hw={1.2:15,2.0:20,6.0:30}
        #  this loop is to run the script on a sample of images and get the 
        #  store the result in a file 
        #
        if 'APERWID' in hdr0:
            if float(hdr0['APERWID'])<=1.3:
                _hw=hw[1.2]
            elif float(hdr0['APERWID'])>1.3 and float(hdr0['APERWID'])<=3.0:
                _hw=hw[2.0]
            elif float(hdr0['APERWID'])>=5:
                _hw=hw[6.0]
            else:
                _hw=20
        else:
            _hw=20
        #####################################################
        _show = False
        # find slit position
        xpos,ypos,aa,bb = floacqastrodef.findcenter(img,_hw,_show)
        pl.draw()
        pl.clf()    

        f[0].header['xpos'] = xpos
        f[0].header['ypos'] = ypos
        f.flush()  # changes are written back to original.fits
        f.close()  # closing the file will also flush an
    else:
        print 'new slit position found'
        ypos = f[0].header['YPOS']
        xpos = f[0].header['XPOS']
        f.close()
    ################
    

    # do astrometry
    os.system('./floacq.py '+img)    
    X0, hdr0 = pyfits.getdata(img, header=True)
    wcs = pywcs.WCS(hdr0)
    # take object position from heder
    objra =  hdr0['CAT-RA']
    objdec = hdr0['CAT-DEC']
    # transform coordinate
    objra, objdec = floacqastrodef.deg2HMS(objra, objdec)
    objcoo = zip([objra],[objdec])
    #  go to pixels
    objpix = wcs.wcs_world2pix(objcoo, 1)
    objxx , objyy = zip(*objpix)
    
    _z1,_z2 = '',''#,lco40.zscale(X0[yy0:yy1,xx0:xx1])
    _z1 = np.percentile(X0,15)
    _z2 = np.percentile(X0,99)
    ax = fig.add_subplot(111)

    if len(X0)<600:
        if hdr0.get('INSTRUME')=='kb37':
            xzero=400
            y0=220
            y1=273
        else:
            xzero=391
            y0=240
            y1=300
        rr=1
    else:
        if hdr0.get('INSTRUME')=='kb37':
            xzero=800
            y0=450
            y1=650
        else:
            xzero=782
            y0=450
            y1=650
        rr=2  # 2
        
    data1 = X0[y0:y1,xzero-50*rr:xzero+50*rr]      #  cut array close to slit position 

    xref0 = xzero-50*rr
    xref1 = xzero+50*rr
    yref0 = y0
    yref1 = y1
    
    ax.imshow(data1, interpolation='nearest', cmap='Greys')#cmap=pl.cm.ocean)
    if xpos:
        ax.plot([xpos-xref0],[ypos-yref0] ,'o', markersize = 10, markeredgecolor = 'r', mfc='none',mew=2)

    ax.plot(spxx-xref0,spyy-yref0 ,'d', markersize = 10, markeredgecolor = 'g', mfc='none',mew=2)
    ax.plot(objxx[0]-xref0,objyy[0]-yref0 ,'s', markersize = 10, markeredgecolor = 'b', mfc='none',mew=2)
    pl.draw()

#    ax.imshow(X0, cmap='gray_r', aspect='equal', interpolation='nearest', origin='lower',vmin=_z1, vmax=_z2)
#    ax.set_xlim(300,500)
#    ax.set_ylim(200,300)
#    ax.plot(spxx,spyy ,'o', markersize = 10, markeredgecolor = 'g', mfc='none',mew=2)
#    ax.plot(objxx,objyy ,'o', markersize = 10, markeredgecolor = 'r', mfc='none',mew=2)
#    if xpos:
#        ax.plot([xpos],[ypos] ,'o', markersize = 10, markeredgecolor = 'b', mfc='none',mew=2)
#    pl.draw()
    fig.savefig(img.replace('.fits', 'slit.png'),dpi= 200)
    pl.clf()

    
os.system('rm output.mp4')

comand = 'cat *slit.png | ffmpeg -f image2pipe -r 5 -i  - output.mp4'
os.system(comand)
#if os.path.isfile('output.mp4'):
#    os.system('rm *slit.png')
    
