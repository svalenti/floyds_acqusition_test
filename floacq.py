#!/usr/bin/env python

import os
import string
import string,sys,re,os
from optparse import OptionParser
from numpy import sqrt
from datetime import datetime
import time
import floacqastrodef
from  floacqastrodef import readkey3, readhdr, updateheader,readlist
description="> Fast astrometry of floyds acq images " 
usage= "%prog  listfile"

def xpa(arg):
    import subprocess
    subproc = subprocess.Popen('xpaset -p ds9 '+arg,shell=True)
    subproc.communicate()

if __name__ == "__main__":
     parser = OptionParser(usage=usage,description=description)
     parser.add_option("-f", "--fit",dest="fitgeo",default='rxyscale',type="str",
                  help='fit fitgeo \t rxyscale, xyscale, shift, general \t[%default]')
     parser.add_option("-s", "--system",dest="system",default='filter',type="str",
                  help='reference photometric system: landolt,sloan,2mass,filter(selected by filter) \t[%default]')
     parser.add_option("-c", "--catalogue",dest="catalogue",default='inst',type="str",
                  help='reference catalogue: usnoa2, usnob2, 2mass,sloan,inst(selected by instrument) \t[%default]')
     parser.add_option("-m", "--method",dest="method",default='vizir',type="str",
                  help='method to query the catalogue: iraf, vizir  \t[%default]')
     parser.add_option("-n", "--number1",dest="number1",default=50,type="int",
                  help='number stars iteraction 1 \t [%default]')
     parser.add_option("-N", "--number2",dest="number2",default=50,type="int",
                  help='number stars iteraction 2 \t [%default]')
     parser.add_option("-x", "--xmatch",dest="match",default=20,type="int",
                  help='minimum number of stars to beleive astronomy \t [%default]')
     parser.add_option("-M", "--number3",dest="number3",default=100,type="int",
                  help='number stars iteraction 3 \t [%default]')
     parser.add_option("-t", "--tollerance1",dest="tollerance1",default=100,type="int",
                  help='tollerance (large) to cross match stars and catalogues  \t [%default]')
     parser.add_option("-T", "--tollerance2",dest="tollerance2",default=30,type="int",
                  help='tollerance (small) to cross match stars and catalogues \t [%default]')
     parser.add_option("-i", "--interactive",dest="interactive",action="store_true")
     parser.add_option("-z", "--zeropoint",dest="zeropoint",help='do zeropoint',action="store_true")
     parser.add_option("-r", "--redo",dest="redo",help='redo astrometry',action="store_true")
     parser.add_option("--xshift",dest="xshift",default=0,type="int",
                  help='x shift in the guess astrometry \t [%default]')
     parser.add_option("--yshift",dest="yshift",default=0,type="int",
                  help='y shift in the guess astrometry \t [%default]')


     option,args = parser.parse_args()
     _method = option.method
     if _method not in ['iraf','vizir']:           sys.argv.append('--help')
     if len(args)<1:                               sys.argv.append('--help')
     option,args = parser.parse_args()
     _start=time.time()
     imglist=readlist(args[0])
     _interactive=option.interactive
     if not _interactive:   _interactive=False 
     _redo=option.redo
     if not _redo:   _redo=False 
     _zeropoint=option.zeropoint
     if not _zeropoint: _zeropoint=False
     _fitgeo = option.fitgeo
     _system = option.system
     _catalogue = option.catalogue
     number1 = option.number1
     number2 = option.number2
     number3 = option.number3
     _t1=option.tollerance1
     _t2=option.tollerance2
     _minnum=option.match
     _xshift=option.xshift
     _yshift=option.yshift
     if _catalogue in ['inst']:  _catalogue2=['2mass','usnoa2']
     #elif _catalogue in ['2mass','usnoa2','usnob1']:   _catalogue2=_catalogue
     else:   _catalogue2=[_catalogue]
     for img in imglist:
       hdr=floacqastrodef.readhdr(img)
       #_wcserr=floacqastrodef.readkey3(hdr,'wcserr')
       _astromet=floacqastrodef.readkey3(hdr,'ASTROMET')
       done=0
       #if float(_wcserr)==0: done=1
       #if done==0:
       if _astromet:
               if float(string.split(_astromet)[0])<2 and float(string.split(_astromet)[1])<2: done=1
       if done==1 and _redo==False:  print '\n#####  astrometry already done'
       else:
         _instrume=floacqastrodef.readkey3(hdr,'instrume')
         if _instrume in ['FLOYDSAG-kb42','FLOYDSAG-kb37','kb42','kb37']:     # FT field are small, half number of star
             number1,number2,number3=number1/2,number2/2,number3/2  
         sexvec=floacqastrodef.sextractor(img)
         try:
            for cat in _catalogue2: #['2mass','usnoa2']:#,'usnob1']:
               rmsx3,rmsy3,num3,fwhmgess,ellgess,ccc,rasys3,decsys3,mbkg3=floacqastrodef.lscastroloop([img],cat,_interactive,number1,number2,number3,'rxyscale',_t1,_t2,sexvec,True,_minnum,_method,_xshift,_yshift)
               if rmsx3 <= 2 and rmsy3 <= 2: break
            print 'no good solution with '+str(_minnum)+' '+str(number1)
            if rmsx3 > 2 and rmsy3 > 2:
               for cat in _catalogue2: #['2mass','usnoa2']:#,'usnob1']:
                   rmsx3,rmsy3,num3,fwhmgess,ellgess,ccc,rasys3,decsys3,mbkg3=floacqastrodef.lscastroloop([img],cat,_interactive,int(number1/2.),int(number2/2.),\
                                                                                                                          int(number3/2.),'rxyscale',_t1,_t2,sexvec,True,int(_minnum/2.),_method,_xshift,_yshift)
                   if rmsx3 <= 2 and rmsy3 <= 2: break
               print 'no good solution with '+str(_minnum/2.)+' '+str(number1/2.)
               if rmsx3 > 2 and rmsy3 > 2:
                   for cat in _catalogue2: #['2mass','usnoa2']:#,'usnob1']:
                       rmsx3,rmsy3,num3,fwhmgess,ellgess,ccc,rasys3,decsys3,mbkg3=floacqastrodef.lscastroloop([img],cat,_interactive,int(10),int(10),\
                                                                                                                          int(25),'rxyscale',_t1,_t2,sexvec,True,int(3),_method,_xshift,_yshift)
            astrostring=str(rmsx3)+' '+str(rmsy3)+' '+str(num3)
            floacqastrodef.updateheader(img,0,{'ASTROMET':[astrostring,'rmsx rmsy nstars']})
         except Exception, e:    
             print e
             rmsx3,rmsy3,num3,fwhmgess,ellgess,ccc,rasys3,decsys3,mbkg3='','','','','','','','',''
             print '\n### problem with astrometry, do you have network ? '
         print img
         if fwhmgess and fwhmgess < 99: 
            print '\n### check astrometry: fine \n### rmsx rmsy nstars: '+astrostring
            floacqastrodef.updateheader(img,0,{'PSF_FWHM':[fwhmgess,'FHWM (arcsec) - computed with sectractor'],'ELLIPTIC':[ellgess,'ellipticity of point sources (1-b/a)'],\
                       'CRDER1':[(1/sqrt(2.))*float(rmsx3),'Random error in axis 1'],'CRDER2':[(1/sqrt(2.))*float(rmsy3),'Random error in axis 2'],\
                       'CUNIT1':['deg','unit of the coord. trans.'],'CUNIT2':['deg','unit of the coord. trans.'],\
                       'CSYER1':[rasys3,'Systematic error (RA_m - Ra_ref)'],'CSYER2':[decsys3,'Systematic error (DEC_m - DEC_ref)']})
       hdr=floacqastrodef.readhdr(img)
       _astromet=floacqastrodef.readkey3(hdr,'ASTROMET')
     _stop=time.time()
     print '\n###time to run '+str((_stop-_start)/60.)
#################################################################################
