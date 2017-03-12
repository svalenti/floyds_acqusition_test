from astropy.io import fits as pyfits
import numpy as np
#########################################################################
def readlist(listfile):
    import string,os,sys,re,glob
    from pyfits import open as opn
    if '*' in listfile:
        imglist=glob.glob(listfile)
    elif ',' in listfile: imglist = string.split(listfile,sep=',')
    else:
        try:            hdulist= opn(listfile)
        except:           hdulist=[]
        if hdulist:            imglist = [listfile]
        else:
           try:
              ff = open(listfile,'r')
              files = ff.readlines()
              ff.close()
              imglist = []
              for ff in files: 
                 ff=re.sub(' ','',ff)
                 if not ff=='\n' and ff[0]!='#':
                    ff=re.sub('\n','',ff)
                    try:
                       hdulist= opn(ff)
                       imglist.append(ff)
                    except:
                       try:
                          correctcard(ff)
                          hdulist= opn(ff)
                          imglist.append(ff)
                       except:                          pass
           except:              sys.exit('\n##### Error ###\n file '+str(listfile)+' do not  exist\n')
    if len(imglist)==0:
           sys.exit('\n##### Error ###\nIf "'+str(listfile)\
                                +'" is an image, it is corrupted \n or is not a list of image\n')
    return imglist
##############################################################################
###############################################################
def readhdr(img):
   from pyfits import open as popen
   try:    hdr=popen(img)[0].header
   except:
      try: 
         correctcard(img)
      except: 
         import sys
         sys.exit('image '+str(img)+' is corrupted, delete it and start again')
      hdr=popen(img)[0].header
   return hdr

#################################################################
################################################
def correctcard(img):
    from  pyfits import open as popen
    from numpy  import asarray
    import re
    hdulist=popen(img)
    a=hdulist[0]._verify('fix')    
    _header=hdulist[0].header
    for i in range(len(a)):
        if not a[i]:
            a[i]=['']
    ww=asarray([i for i in range(len(a)) if (re.sub(' ','',a[i][0])!='')])
    if len(ww)>0:
        newheader=[]
        headername=[]
        for j in _header.items():
            headername.append(j[0])
            newheader.append(j[1])
        hdulist.close()
        imm=popen(img,mode='update')
        _header=imm[0].header
        for i in ww:
            if headername[i]:
                try:
                    _header.update(headername[i],newheader[i])
                except:
                    _header.update(headername[i],'xxxx')
        imm.flush()
        imm.close()
######################################################################################################


def readkey3(hdr,keyword):
    import re,string,sys
    import pyfits
    if pyfits.__version__:
         if int(re.sub('\.','',str(pyfits.__version__))[:2])<=30:  aa='HIERARCH '
         else: aa=''
    else:  aa=''
    try:    _instrume=hdr.get('INSTRUME').lower()
    except: _instrume='none'
    if _instrume in ['FLOYDSAG-kb42','FLOYDSAG-kb37','kb37']:
         if not hdr.get('HDRVER'):
              useful_keys = {'object'    : 'OBJECT',\
                             'date-obs'  : 'DATE-OBS',\
                             'ut'        : 'UTSTART',\
                             'obstype'   : 'OBSTYPE',\
                             'RA'        : 'RA',\
                             'DEC'       : 'DEC',\
                             'datamin'   : -100,\
                             'datamax'   : 60000,\
                             'grpid'     : 'GRPUID',\
                             'exptime'   : 'EXPTIME',\
                             'JD'        : 'MJD',\
                             'lamp'      : 'LMP_ID',\
                             'gain'      : 'GAIN',\
                             'instrume'  : 'INSTRUME',\
                             'grism'     : 'GRISM',\
                             'ron'       : 'RDNOISE',\
                             'airmass'   : 'AIRMASS',\
                             'slit'      : 'APERWID',\
                             'telescop'  : 'TELESCOP'}
         else:
              useful_keys = {'object'    : 'OBJECT',\
                             'date-obs'  : 'DATE-OBS',\
                             'ut'        : 'UTSTART',\
                             'obstype'   : 'OBSTYPE',\
                             'RA'        : 'RA',\
                             'DEC'       : 'DEC',\
                             'datamin'   : -100,\
                             'datamax'   : 60000,\
                             'grpid'     : 'GRPUID',\
                             'exptime'   : 'EXPTIME',\
                             'JD'        : 'MJD-OBS',\
                             'lamp'      : 'LMP1ID',\
                             'gain'      : 'GAIN',\
                             'instrume'  : 'INSTRUME',\
                             'grism'     : 'GRISM',\
                             'ron'       : 'RDNOISE',\
                             'airmass'   : 'AIRMASS',\
                             'slit'      : 'APERWID',\
                             'telescop'  : 'TELESCOP'}              
    else: 
          useful_keys = {'object'    : 'OBJECT',\
                         'date-obs'  : 'DATE-OBS'}
    if keyword in useful_keys:
        if type(useful_keys[keyword])==float:
            value=useful_keys[keyword]
        else:
            value=hdr.get(useful_keys[keyword])
            if keyword=='date-obs':
                import string,re
                try:
                   value=re.sub('-','',string.split(value,'T')[0])
                except:
                   pass
            elif keyword=='ut':
                import string,re
                try:
                   value=string.split(value,'T')[1]
                except:
                   pass
            elif keyword=='JD':       value=float(value)+0.5
            elif keyword=='instrume':      value=value.lower()
            elif keyword=='grism':
                 if not value: value='full'
            elif keyword=='RA':
                 import string,re
                 value0=string.split(value,':')
                 value=((float(value0[0])+((float(value0[1])+(float(value0[2])/60.))/60.))*15)
            elif keyword=='DEC':
                 import string,re
                 value0=string.split(value,':')
                 if '-' in str(value0[0]):
                      value=((-1)*(abs(float(value0[0]))+((float(value0[1])+(float(value0[2])/60.))/60.)))
                 else:
                      value=(float(value0[0])+((float(value0[1])+(float(value0[2])/60.))/60.))
            elif keyword=='slit':     
                 value=re.sub('\"','',re.sub('slit','',str(value)))
            elif keyword=='object':
                 value=re.sub('\}','',value)
                 value=re.sub('\{','',value)
                 value=re.sub('\[','',value)
                 value=re.sub('\]','',value)
                 value=re.sub('\(','',value)
                 value=re.sub('\)','',value)
                 value=re.sub('-','',value)
                 value=re.sub(' ','',value)
    else:
       if keyword=='date-night':
            import datetime
            _date=readkey3(hdr,'DATE-OBS')
            a=(datetime.datetime.strptime(string.split(_date,'.')[0],"20%y-%m-%dT%H:%M:%S")-datetime.timedelta(.0)).isoformat()
            value=re.sub('-','',string.split(a,'T')[0])
       elif keyword=='TELID':
            value=hdr.get(keyword)
            value=re.sub('-','',value)
            if value in ['fts','2m0b']: value='fts'
            elif value in ['ftn','2m0a']: value='ftn'
            else:  sys.exit('Warning: keyword not valid')
       else:
          try:     value=hdr.get(keyword)
          except:       sys.exit('Warning: keyword not valid')
    if type(value) == str:    value=re.sub('\#','',value)
    return value

#######################################################

######################################################################################################
def updateheader(image,dimension,headerdict):
    from pyfits import open as opp
    try:
        imm=opp(image,mode='update')
        _header=imm[dimension].header
################################
#   change way to update to speed up the process
#   now get dictionary   08 12  2012
################################
        for i in headerdict.keys():
           _header.update(i,headerdict[i][0],headerdict[i][1])
###################################################
#        _header.update(_headername,_value,commento)
        imm.flush()
        imm.close()
    except:
        print 'warning: problem to update header, try to correct header format ....'
        correctcard(image)
        try:
            imm=opp(image,mode='update')
            _header=imm[dimension].header
###################################################
            for i in headerdict.keys():
               _header.update(i,headerdict[i][0],headerdict[i][1])
###################################################
            _header.update(_headername,_value,commento)
            imm.flush()
            imm.close()
        except:
           print 'niente'
#           import sys
#            sys.exit('error: not possible update header')
#################################################################################################


def wcsstart(img,CRPIX1='',CRPIX2=''):
    from numpy import pi, sin, cos 
    import floacqastrodef
    from floacqastrodef import readhdr,readkey3
    hdr=floacqastrodef.readhdr(img)
    _instrume=readkey3(hdr,'instrume')
    _RA=floacqastrodef.readkey3(hdr,'RA')
    _DEC=floacqastrodef.readkey3(hdr,'DEC')
    _RA,_DEC=floacqastrodef.deg2HMS(_RA,_DEC)
    _xdimen=floacqastrodef.readkey3(hdr,'NAXIS1')
    _ydimen=floacqastrodef.readkey3(hdr,'NAXIS2')
    _CCDXBIN=floacqastrodef.readkey3(hdr,'CCDXBIN')
    if _instrume in ['FLOYDSAG-kb42','FLOYDSAG-kb37','kb37']:
        angle=83.4-floacqastrodef.readkey3(hdr,'ROTSKYPA')#posang)
        theta=(angle*pi/180.)
        pixscale=0.253   # 0.253
        CDELT0=pixscale/3600.
        CD1_1=(-1)*CDELT0*cos(theta)
        CD1_2=(-1)*CDELT0*sin(theta)
        CD2_1=(-1)*CDELT0*sin(theta)
        CD2_2=CDELT0*cos(theta)
        if not CRPIX1:        CRPIX1= 750.
        else: CRPIX1= 750.+CRPIX1
        if not CRPIX2:        CRPIX2= 600.
        else: CRPIX2= 600.+CRPIX2
        CDELT1=2
        CDELT2=2
    elif _instrume in ['kb74','kb76']:
        angle=readkey3(hdr,'ROLLERDR')#posang)
        theta=(angle*pi/180.)
        CDELT0=0.000129722   # 1.3042840792028E-4   8.43604528922325E-5  #6.6888889999999995e-05
        CD1_1=(-1)*CDELT0*cos(theta)
        CD2_2=CDELT0*cos(theta)
        CD1_2=abs(CDELT0)*(abs(CDELT0)/CDELT0)*sin(theta)
        CD2_1=abs(CDELT0)*(abs(CDELT0)/CDELT0)*sin(theta)
        if not CRPIX1:        CRPIX1= readkey3(hdr,'ROTCENTX')
        else: CRPIX1= 1000.+CRPIX1
        if not CRPIX2:        CRPIX2= readkey3(hdr,'ROTCENTY')
        else: CRPIX2= 1000.+CRPIX2
        CDELT1=2
        CDELT2=2
    elif _instrume in ['fs01','fs02']:
        angle=readkey3(hdr,'ROTSKYPA')#posang)
        theta=(angle*pi/180.)
 #       pixscale=0.30*_CCDXBIN
 #       CDELT0=pixscale/3600.
        CDELT0=0.000083568667  # 8.43604528922325E-5  #6.6888889999999995e-05
        CD1_1=(-1)*CDELT0*cos(theta)
        CD2_2=CDELT0*cos(theta)
        CD1_2=abs(CDELT0)*(abs(CDELT0)/CDELT0)*sin(theta)
        CD2_1=abs(CDELT0)*(abs(CDELT0)/CDELT0)*sin(theta)
        if not CRPIX1:        CRPIX1= 1024.
        else: CRPIX1= 1024.+CRPIX1
        if not CRPIX2:        CRPIX2= 1024.
        else: CRPIX2= 1024.+CRPIX2
        CDELT1=0.000083568694*(-1)
        CDELT2=0.000083568694
    elif _instrume in ['fs03']:
        angle=readkey3(hdr,'ROTSKYPA')#posang)
        theta=(angle*pi/180.)
 #       pixscale=0.30*_CCDXBIN
 #       CDELT0=pixscale/3600.
        CDELT0=0.000083705976   #8.43604528922325E-5  #6.6888889999999995e-05
        CD1_1=(-1)*CDELT0*cos(theta)
        CD2_2=CDELT0*cos(theta)
        CD1_2=abs(CDELT0)*(abs(CDELT0)/CDELT0)*sin(theta)
        CD2_1=abs(CDELT0)*(abs(CDELT0)/CDELT0)*sin(theta)
        if not CRPIX1:        CRPIX1= 1024.
        else: CRPIX1= 1024.+CRPIX1
        if not CRPIX2:        CRPIX2= 1024.
        else: CRPIX2= 1024.+CRPIX2
        CDELT1=0.000083705976*(-1)
        CDELT2=0.000083705976
    elif _instrume.lower() in ['em03','em01']:
        theta=(angle*pi/180.)
        CDELT0=7.63077724258886e-05 #7.7361111111111123e-05 
        CD1_1=(-1)*CDELT0*cos(theta)
        CD2_2=CDELT0*cos(theta)
        CD1_2=abs(CDELT0)*(abs(CDELT0)/CDELT0)*sin(theta)
        CD2_1=abs(CDELT0)*(abs(CDELT0)/CDELT0)*sin(theta)
        if not CRPIX1:        CRPIX1= readkey3(hdr,'ROTCENTX')
        if not CRPIX2:        CRPIX2= readkey3(hdr,'ROTCENTY')-100
        CDELT1=2
        CDELT2=2
    else:  print '\n### ERROR: instument not found !!!'

    CTYPE1  = 'RA---TAN'  
    CTYPE2  = 'DEC--TAN' 
    CRVAL1=_RA
    CRVAL2=_DEC
    WCSDIM  =                   2  
    LTM1_1  =                   1. 
    LTM2_2  =                   1.
    WAT0_001= 'system=image'
    WAT1_001= 'wtype=tan axtype=ra'
    WAT2_001= 'wtype=tan axtype=dec'
    updateheader(img,0,{'CTYPE1':[CTYPE1,''], 'CTYPE2':[CTYPE2,''], 'CRVAL1':[CRVAL1,''], 'CRVAL2':[CRVAL2,''],\
                'CRPIX1':[CRPIX1,''], 'CRPIX2':[CRPIX2,''], 'CDELT1':[CDELT1,''], 'CDELT2':[CDELT2,''],\
                'CD1_1':[CD1_1,''], 'CD2_2':[CD2_2,''], 'CD1_2':[CD1_2,''], 'CD2_1':[CD2_1,''],\
                'WCSDIM':[WCSDIM,'']})

##########################################################

def sextractor(img):
        import os
        from pyraf import iraf
        from iraf import proto
        from numpy import compress,array,asarray
        
#        namesex=defsex('default.sex')
        os.system('sex '+img+' -c '+namesex+' > _logsex')
        delete(namesex)
        delete('_logsex')
        xpix=iraf.proto.fields('detections.cat',fields='2',Stdout=1)
        ypix=iraf.proto.fields('detections.cat',fields='3',Stdout=1)
        cm=iraf.proto.fields('detections.cat',fields='4',Stdout=1)
        cl=iraf.proto.fields('detections.cat',fields='7',Stdout=1)
        fw=iraf.proto.fields('detections.cat',fields='8',Stdout=1)
        ell=iraf.proto.fields('detections.cat',fields='9',Stdout=1)
        bkg=iraf.proto.fields('detections.cat',fields='10',Stdout=1)

        cl=compress((array(xpix)!=''),array(cl,float))
        cm=compress((array(xpix)!=''),array(cm,float))
        fw=compress((array(xpix)!=''),array(fw,float))
        ell=compress((array(xpix)!=''),array(ell,float))
        bkg=compress((array(xpix)!=''),array(bkg,float))
        ypix=compress((array(xpix)!=''),array(ypix,float))
        xpix=compress((array(xpix)!=''),array(xpix,float))
        try:
            ww=asarray([i for i in range(len(xpix)) if ((xpix[i]<2010) or (ypix[i]<2010))])
            cl,cm,fw,ell,xpix,ypix,bkg=cl[ww],cm[ww],fw[ww],ell[ww],xpix[ww],ypix[ww],bkg[ww]

            ww=asarray([i for i in range(len(xpix)) if ((xpix[i]>20) or (ypix[i]<2010))])
            cl,cm,fw,ell,xpix,ypix,bkg=cl[ww],cm[ww],fw[ww],ell[ww],xpix[ww],ypix[ww],bkg[ww]

            ww=asarray([i for i in range(len(xpix)) if (xpix[i]>3)])
            cl,cm,fw,ell,xpix,ypix,bkg=cl[ww],cm[ww],fw[ww],ell[ww],xpix[ww],ypix[ww],bkg[ww]

            cl=compress((array(fw)<=15)&(array(fw)>=-2),array(cl))
            cm=compress((array(fw)<=15)&(array(fw)>=-2),array(cm))
            xpix=compress((array(fw)<=15)&(array(fw)>=-2),array(xpix))
            ypix=compress((array(fw)<=15)&(array(fw)>=-2),array(ypix))
            ell=compress((array(fw)<=15)&(array(fw)>=-2),array(ell))
            bkg=compress((array(fw)<=15)&(array(fw)>=-2),array(bkg))
            fw=compress((array(fw)<=15)&(array(fw)>=-2),array(fw))
        except: 
            xpix,ypix,fw,cl,cm,ell=[],[],[],[],[],[]
            print '\n### ERROR Filtering the sextractor detections, please check that sextractor is working ......'
        delete('detections.cat')
        return xpix,ypix,fw,cl,cm,ell,bkg

###############################################################################

def defsex(namefile):
    import lsc
    import string,re,os
    sexfile=lsc.__path__[0]+'/standard/sex/default.sex'
    f=open(sexfile,'r')
    ss=f.readlines()
    f.close()   
    ff=open(namefile,'w')
    for i in ss:
        if string.count(i,'PARAMETERS_NAME')==1:
            ff.write('PARAMETERS_NAME  "'+lsc.__path__[0]+'/standard/sex/default.param"\n')
        elif string.count(i,'FILTER_NAME')==1:
            ff.write('FILTER_NAME  "'+lsc.__path__[0]+'/standard/sex/default.conv"\n')
        elif string.count(i,'STARNNW_NAME')==1:
            ff.write('STARNNW_NAME "'+lsc.__path__[0]+'/standard/sex/default.nnw"\n')
        else:
            ff.write(i)
    ff.close()
    return namefile

###############################################################################################

def vizq(_ra,_dec,catalogue,radius):
    ''' Query vizquery '''
    import os,string,re
    import floacqastrodef
    from numpy import array,compress
#    _site='vizier.cfa.harvard.edu'
    _site='vizier.u-strasbg.fr'
    cat={'usnoa2':['I/252/out','USNO-A2.0','USNO-A2.0,Rmag'],\
         '2mass':['II/246/out','2MASS','2MASS,Jmag'],\
         'landolt':['II/183A/table2','','Vmag,B-V,U-B,V-R,R-I,Star,e_Vmag'],\
         'apass':['I/322A/out','','Bmag,Vmag,gmag,rmag,imag,e_Vmag,e_Bmag,e_gmag,e_rmag,e_imag,UCAC4'],\
         'usnob1':['I/284/out','USNO-B1.0','R2mag'],'sdss7':['II/294/sdss7','','objID,umag,gmag,rmag,imag,zmag,gc'],\
         'sdss9':['V/139/sdss9','','objID,umag,gmag,rmag,imag,zmag,e_umag,e_gmag,e_rmag,e_imag,e_zmag,gc'],\
         'sdss7':['II/294/sdss7','','objID,umag,gmag,rmag,imag,zmag,e_umag,e_gmag,e_rmag,e_imag,e_zmag,gc'],\
         'sdss8':['II/306/sdss8','','objID,umag,gmag,rmag,imag,zmag,e_umag,e_gmag,e_rmag,e_imag,e_zmag,gc']}

    a=os.popen('vizquery -mime=tsv  -site='+_site+' -source='+cat[catalogue][0]+\
                   ' -c.ra='+str(_ra)+' -c.dec='+str(_dec)+' -c.eq=J2000 -c.rm='+str(radius)+\
                   ' -c.geom=b -oc.form=h -sort=_RA*-c.eq -out.add=_RAJ2000,_DEJ2000 -out.max=10000 -out='+\
                   cat[catalogue][1]+' -out='+cat[catalogue][2]+'').read()
    print 'vizquery -mime=tsv  -site='+_site+' -source='+cat[catalogue][0]+\
                   ' -c.ra='+str(_ra)+' -c.dec='+str(_dec)+' -c.eq=J2000 -c.rm='+str(radius)+\
                   ' -c.geom=b -oc.form=h -sort=_RA*-c.eq -out.add=_RAJ2000,_DEJ2000 -out.max=10000 -out='+\
                   cat[catalogue][1]+' -out='+cat[catalogue][2]+''
    aa=a.split('\n')
    bb=[]
    for i in aa:
        if i and i[0]!='#':   bb.append(i)
    _ra,_dec,_name,_mag=[],[],[],[]
    for ii in bb[3:]:
        aa=ii.split('\t')
        rr,dd=floacqastrodef.deg2HMS(ra=re.sub(' ',':',aa[0]), dec=re.sub(' ',':',aa[1]), round=False)
        _ra.append(rr)
        _dec.append(dd)
        _name.append(aa[2])
    dictionary={'ra':_ra,'dec':_dec,'id':_name}
    sss=string.split(cat[catalogue][2],',')
    for ii in sss: dictionary[ii]=[]
    for ii in bb[3:]:
        aa=ii.split('\t')
        for gg in range(0,len(sss)):
           if sss[gg] not in ['UCAC4','id']:
              try:
                 dictionary[sss[gg]].append(float(aa[2+gg]))
              except:    
                 dictionary[sss[gg]].append(float(9999))
           else:
                 dictionary[sss[gg]].append(str(aa[2+gg]))

    if catalogue in ['sdss7','sdss9','sdss8']:
        dictionary['u']=dictionary['umag']
        dictionary['g']=dictionary['gmag']
        dictionary['r']=dictionary['rmag']
        dictionary['i']=dictionary['imag']
        dictionary['z']=dictionary['zmag']
        dictionary['uerr']=dictionary['e_umag']
        dictionary['gerr']=dictionary['e_gmag']
        dictionary['rerr']=dictionary['e_rmag']
        dictionary['ierr']=dictionary['e_imag']
        dictionary['zerr']=dictionary['e_zmag']
        for key in dictionary.keys():
           if key!='r':
              dictionary[key]=compress((array(dictionary['r'])<19)&(array(dictionary['r']>13)),dictionary[key])
        dictionary['r']=compress((array(dictionary['r'])<19)&(array(dictionary['r']>13)),dictionary['r'])

    elif  catalogue=='landolt':
        dictionary['B']=array(dictionary['Vmag'])+array(dictionary['B-V'])
        dictionary['U']=array(dictionary['B'])+array(dictionary['U-B'])
        dictionary['V']=array(dictionary['Vmag'])
        dictionary['Verr']=array(dictionary['e_Vmag'])
        dictionary['R']=array(dictionary['Vmag'])-array(dictionary['V-R'])
        dictionary['I']=array(dictionary['R'])-array(dictionary['R-I'])
        dictionary['id']=array(dictionary['Star'])
    elif  catalogue=='apass':
        dictionary['B']=array(dictionary['Bmag'])
        dictionary['V']=array(dictionary['Vmag'])
        dictionary['g']=array(dictionary['gmag'])
        dictionary['r']=array(dictionary['rmag'])
        dictionary['i']=array(dictionary['imag'])
        dictionary['Berr']=array(dictionary['e_Bmag'],float)/100.
        dictionary['Verr']=array(dictionary['e_Vmag'],float)/100.
        dictionary['gerr']=array(dictionary['e_gmag'],float)/100.
        dictionary['rerr']=array(dictionary['e_rmag'],float)/100.
        dictionary['ierr']=array(dictionary['e_imag'],float)/100.
        dictionary['id']=array(dictionary['UCAC4'],str)
        for key in dictionary.keys():
           if key!='r':
              dictionary[key]=compress((array(dictionary['r'])<19)&(array(dictionary['r']>10.5)),dictionary[key])
        dictionary['r']=compress((array(dictionary['r'])<19)&(array(dictionary['r']>10.5)),dictionary['r'])
    return dictionary
####################################################################################################

##############################################################################
def delete(listfile):
    import os,string,re,glob
    if listfile[0]=='@':   
        ff = open(listfile[1:])
        files = ff.readlines()
        imglist = []
        for ff in files: 
            ff=re.sub(' ','',ff)
            if not ff=='\n' and ff[0]!='#':
                ff=re.sub('\n','',ff)
                imglist.append(ff)
    elif ',' in listfile: imglist = string.split(listfile,sep=',')
    else:       imglist=[listfile]    
    lista=[]
    for _file in imglist:   lista=lista+glob.glob(_file)
    if lista:
        for _file in lista:
            try:          os.system('rm '+_file)
            except:       pass
###############################################################
#############################################

def deg2HMS(ra='', dec='', round=False):
      import string
      RA, DEC= '', ''
      if dec:
          if string.count(str(dec),':')==2:
              dec00=string.split(dec,':')
              dec0,dec1,dec2=float(dec00[0]),float(dec00[1]),float(dec00[2])
              if '-' in str(dec0):       DEC=(-1)*((dec2/60.+dec1)/60.+((-1)*dec0))
              else:                      DEC=(dec2/60.+dec1)/60.+dec0
          else:
              if str(dec)[0]=='-':      dec0=(-1)*abs(int(dec))
              else:                     dec0=abs(int(dec))
              dec1=int((abs(dec)-abs(dec0))*(60))
              dec2=((((abs(dec))-abs(dec0))*60)-abs(dec1))*60
              DEC='00'[len(str(dec0)):]+str(dec0)+':'+'00'[len(str(dec1)):]+str(dec1)+':'+'00'[len(str(int(dec2))):]+str(dec2)
      if ra:
          if string.count(str(ra),':')==2:
              ra00=string.split(ra,':')
              ra0,ra1,ra2=float(ra00[0]),float(ra00[1]),float(ra00[2])
              RA=((ra2/60.+ra1)/60.+ra0)*15.
          else:
              ra0=int(ra/15.)
              ra1=int(((ra/15.)-ra0)*(60))
              ra2=((((ra/15.)-ra0)*60)-ra1)*60
              RA='00'[len(str(ra0)):]+str(ra0)+':'+'00'[len(str(ra1)):]+str(ra1)+':'+'00'[len(str(int(ra2))):]+str(ra2)
      if ra and dec:          return RA, DEC
      else:                   return RA or DEC

###########################################################################

def lscastroloop(imglist,catalogue,_interactive,number1,number2,number3,_fitgeo,_tollerance1,_tollerance2,sexvec='',_guess=False,_numin=4,method='vizir',xshift=0,yshift=0):
    import floacqastrodef
    from floacqastrodef import lscastrometry2
    from numpy import median, array
    import math
    import datetime
    import time
    _imex=False
    for img in imglist:
        hdr=floacqastrodef.readhdr(img)
        _instrume=floacqastrodef.readkey3(hdr,'instrume')
    if catalogue=='inst' and _instrume in ['FLOYDSAG-kb42','FLOYDSAG-kb37','kb37']:
        angle=readkey3(hdr,'ROLLERDR')#posang)
        theta=(angle*pi/180.)
        CDELT0=0.000129722   # 1.3042840792028E-4   8.43604528922325E-5  #6.6888889999999995e-05
        CD1_1=(-1)*CDELT0*cos(theta)
        CD2_2=CDELT0*cos(theta)
        CD1_2=abs(CDELT0)*(abs(CDELT0)/CDELT0)*sin(theta)
        CD2_1=abs(CDELT0)*(abs(CDELT0)/CDELT0)*sin(theta)
        if not CRPIX1:        CRPIX1= readkey3(hdr,'ROTCENTX')
        else: CRPIX1= 1000.+CRPIX1
        if not CRPIX2:        CRPIX2= readkey3(hdr,'ROTCENTY')
        else: CRPIX2= 1000.+CRPIX2
        CDELT1=2
        CDELT2=2
    elif _instrume in ['fs01','fs02']:
        angle=readkey3(hdr,'ROTSKYPA')#posang)
        theta=(angle*pi/180.)
 #       pixscale=0.30*_CCDXBIN
 #       CDELT0=pixscale/3600.
        CDELT0=0.000083568667  # 8.43604528922325E-5  #6.6888889999999995e-05
        CD1_1=(-1)*CDELT0*cos(theta)
        CD2_2=CDELT0*cos(theta)
        CD1_2=abs(CDELT0)*(abs(CDELT0)/CDELT0)*sin(theta)
        CD2_1=abs(CDELT0)*(abs(CDELT0)/CDELT0)*sin(theta)
        if not CRPIX1:        CRPIX1= 1024.
        else: CRPIX1= 1024.+CRPIX1
        if not CRPIX2:        CRPIX2= 1024.
        else: CRPIX2= 1024.+CRPIX2
        CDELT1=0.000083568694*(-1)
        CDELT2=0.000083568694
    elif _instrume in ['fs03']:
        angle=readkey3(hdr,'ROTSKYPA')#posang)
        theta=(angle*pi/180.)
 #       pixscale=0.30*_CCDXBIN
 #       CDELT0=pixscale/3600.
        CDELT0=0.000083705976   #8.43604528922325E-5  #6.6888889999999995e-05
        CD1_1=(-1)*CDELT0*cos(theta)
        CD2_2=CDELT0*cos(theta)
        CD1_2=abs(CDELT0)*(abs(CDELT0)/CDELT0)*sin(theta)
        CD2_1=abs(CDELT0)*(abs(CDELT0)/CDELT0)*sin(theta)
        if not CRPIX1:        CRPIX1= 1024.
        else: CRPIX1= 1024.+CRPIX1
        if not CRPIX2:        CRPIX2= 1024.
        else: CRPIX2= 1024.+CRPIX2
        CDELT1=0.000083705976*(-1)
        CDELT2=0.000083705976
    elif _instrume.lower() in ['em03','em01']:
        theta=(angle*pi/180.)
        CDELT0=7.63077724258886e-05 #7.7361111111111123e-05 
        CD1_1=(-1)*CDELT0*cos(theta)
        CD2_2=CDELT0*cos(theta)
        CD1_2=abs(CDELT0)*(abs(CDELT0)/CDELT0)*sin(theta)
        CD2_1=abs(CDELT0)*(abs(CDELT0)/CDELT0)*sin(theta)
        if not CRPIX1:        CRPIX1= readkey3(hdr,'ROTCENTX')
        if not CRPIX2:        CRPIX2= readkey3(hdr,'ROTCENTY')-100
        CDELT1=2
        CDELT2=2
    else:  print '\n### ERROR: instument not found !!!'

    CTYPE1  = 'RA---TAN'  
    CTYPE2  = 'DEC--TAN' 
    CRVAL1=_RA
    CRVAL2=_DEC
    WCSDIM  =                   2  
    LTM1_1  =                   1. 
    LTM2_2  =                   1.
    WAT0_001= 'system=image'
    WAT1_001= 'wtype=tan axtype=ra'
    WAT2_001= 'wtype=tan axtype=dec'
    updateheader(img,0,{'CTYPE1':[CTYPE1,''], 'CTYPE2':[CTYPE2,''], 'CRVAL1':[CRVAL1,''], 'CRVAL2':[CRVAL2,''],\
                'CRPIX1':[CRPIX1,''], 'CRPIX2':[CRPIX2,''], 'CDELT1':[CDELT1,''], 'CDELT2':[CDELT2,''],\
                'CD1_1':[CD1_1,''], 'CD2_2':[CD2_2,''], 'CD1_2':[CD1_2,''], 'CD2_1':[CD2_1,''],\
                'WCSDIM':[WCSDIM,'']})




###############################################################################


def vizq(_ra,_dec,catalogue,radius):
    ''' Query vizquery '''
    import os,string,re
    import floacqastrodef
    from numpy import array,compress
#    _site='vizier.cfa.harvard.edu'
    _site='vizier.u-strasbg.fr'
    cat={'usnoa2':['I/252/out','USNO-A2.0','USNO-A2.0,Rmag'],\
         '2mass':['II/246/out','2MASS','2MASS,Jmag'],\
         'landolt':['II/183A/table2','','Vmag,B-V,U-B,V-R,R-I,Star,e_Vmag'],\
         'apass':['I/322A/out','','Bmag,Vmag,gmag,rmag,imag,e_Vmag,e_Bmag,e_gmag,e_rmag,e_imag,UCAC4'],\
         'usnob1':['I/284/out','USNO-B1.0','R2mag'],'sdss7':['II/294/sdss7','','objID,umag,gmag,rmag,imag,zmag,gc'],\
         'sdss9':['V/139/sdss9','','objID,umag,gmag,rmag,imag,zmag,e_umag,e_gmag,e_rmag,e_imag,e_zmag,gc'],\
         'sdss7':['II/294/sdss7','','objID,umag,gmag,rmag,imag,zmag,e_umag,e_gmag,e_rmag,e_imag,e_zmag,gc'],\
         'sdss8':['II/306/sdss8','','objID,umag,gmag,rmag,imag,zmag,e_umag,e_gmag,e_rmag,e_imag,e_zmag,gc']}

    a=os.popen('vizquery -mime=tsv  -site='+_site+' -source='+cat[catalogue][0]+\
                   ' -c.ra='+str(_ra)+' -c.dec='+str(_dec)+' -c.eq=J2000 -c.rm='+str(radius)+\
                   ' -c.geom=b -oc.form=h -sort=_RA*-c.eq -out.add=_RAJ2000,_DEJ2000 -out.max=10000 -out='+\
                   cat[catalogue][1]+' -out='+cat[catalogue][2]+'').read()
    print 'vizquery -mime=tsv  -site='+_site+' -source='+cat[catalogue][0]+\
                   ' -c.ra='+str(_ra)+' -c.dec='+str(_dec)+' -c.eq=J2000 -c.rm='+str(radius)+\
                   ' -c.geom=b -oc.form=h -sort=_RA*-c.eq -out.add=_RAJ2000,_DEJ2000 -out.max=10000 -out='+\
                   cat[catalogue][1]+' -out='+cat[catalogue][2]+''
    aa=a.split('\n')
    bb=[]
    for i in aa:
        if i and i[0]!='#':   bb.append(i)
    _ra,_dec,_name,_mag=[],[],[],[]
    for ii in bb[3:]:
        aa=ii.split('\t')
        rr,dd=floacqastrodef.deg2HMS(ra=re.sub(' ',':',aa[0]), dec=re.sub(' ',':',aa[1]), round=False)
        _ra.append(rr)
        _dec.append(dd)
        _name.append(aa[2])
    dictionary={'ra':_ra,'dec':_dec,'id':_name}
    sss=string.split(cat[catalogue][2],',')
    for ii in sss: dictionary[ii]=[]
    for ii in bb[3:]:
        aa=ii.split('\t')
        for gg in range(0,len(sss)):
           if sss[gg] not in ['UCAC4','id']:
              try:
                 dictionary[sss[gg]].append(float(aa[2+gg]))
              except:    
                 dictionary[sss[gg]].append(float(9999))
           else:
                 dictionary[sss[gg]].append(str(aa[2+gg]))

    if catalogue in ['sdss7','sdss9','sdss8']:
        dictionary['u']=dictionary['umag']
        dictionary['g']=dictionary['gmag']
        dictionary['r']=dictionary['rmag']
        dictionary['i']=dictionary['imag']
        dictionary['z']=dictionary['zmag']
        dictionary['uerr']=dictionary['e_umag']
        dictionary['gerr']=dictionary['e_gmag']
        dictionary['rerr']=dictionary['e_rmag']
        dictionary['ierr']=dictionary['e_imag']
        dictionary['zerr']=dictionary['e_zmag']
        for key in dictionary.keys():
           if key!='r':
              dictionary[key]=compress((array(dictionary['r'])<19)&(array(dictionary['r']>13)),dictionary[key])
        dictionary['r']=compress((array(dictionary['r'])<19)&(array(dictionary['r']>13)),dictionary['r'])

    elif  catalogue=='landolt':
        dictionary['B']=array(dictionary['Vmag'])+array(dictionary['B-V'])
        dictionary['U']=array(dictionary['B'])+array(dictionary['U-B'])
        dictionary['V']=array(dictionary['Vmag'])
        dictionary['Verr']=array(dictionary['e_Vmag'])
        dictionary['R']=array(dictionary['Vmag'])-array(dictionary['V-R'])
        dictionary['I']=array(dictionary['R'])-array(dictionary['R-I'])
        dictionary['id']=array(dictionary['Star'])
    elif  catalogue=='apass':
        dictionary['B']=array(dictionary['Bmag'])
        dictionary['V']=array(dictionary['Vmag'])
        dictionary['g']=array(dictionary['gmag'])
        dictionary['r']=array(dictionary['rmag'])
        dictionary['i']=array(dictionary['imag'])
        dictionary['Berr']=array(dictionary['e_Bmag'],float)/100.
        dictionary['Verr']=array(dictionary['e_Vmag'],float)/100.
        dictionary['gerr']=array(dictionary['e_gmag'],float)/100.
        dictionary['rerr']=array(dictionary['e_rmag'],float)/100.
        dictionary['ierr']=array(dictionary['e_imag'],float)/100.
        dictionary['id']=array(dictionary['UCAC4'],str)
        for key in dictionary.keys():
           if key!='r':
              dictionary[key]=compress((array(dictionary['r'])<19)&(array(dictionary['r']>10.5)),dictionary[key])
        dictionary['r']=compress((array(dictionary['r'])<19)&(array(dictionary['r']>10.5)),dictionary['r'])
    return dictionary
####################################################################################################

##############################################################################
def delete(listfile):
    import os,string,re,glob
    if listfile[0]=='@':   
        ff = open(listfile[1:])
        files = ff.readlines()
        imglist = []
        for ff in files: 
            ff=re.sub(' ','',ff)
            if not ff=='\n' and ff[0]!='#':
                ff=re.sub('\n','',ff)
                imglist.append(ff)
    elif ',' in listfile: imglist = string.split(listfile,sep=',')
    else:       imglist=[listfile]    
    lista=[]
    for _file in imglist:   lista=lista+glob.glob(_file)
    if lista:
        for _file in lista:
            try:          os.system('rm '+_file)
            except:       pass
###############################################################
#############################################

def deg2HMS(ra='', dec='', round=False):
      import string
      RA, DEC= '', ''
      if dec:
          if string.count(str(dec),':')==2:
              dec00=string.split(dec,':')
              dec0,dec1,dec2=float(dec00[0]),float(dec00[1]),float(dec00[2])
              if '-' in str(dec0):       DEC=(-1)*((dec2/60.+dec1)/60.+((-1)*dec0))
              else:                      DEC=(dec2/60.+dec1)/60.+dec0
          else:
              if str(dec)[0]=='-':      dec0=(-1)*abs(int(dec))
              else:                     dec0=abs(int(dec))
              dec1=int((abs(dec)-abs(dec0))*(60))
              dec2=((((abs(dec))-abs(dec0))*60)-abs(dec1))*60
              DEC='00'[len(str(dec0)):]+str(dec0)+':'+'00'[len(str(dec1)):]+str(dec1)+':'+'00'[len(str(int(dec2))):]+str(dec2)
      if ra:
          if string.count(str(ra),':')==2:
              ra00=string.split(ra,':')
              ra0,ra1,ra2=float(ra00[0]),float(ra00[1]),float(ra00[2])
              RA=((ra2/60.+ra1)/60.+ra0)*15.
          else:
              ra0=int(ra/15.)
              ra1=int(((ra/15.)-ra0)*(60))
              ra2=((((ra/15.)-ra0)*60)-ra1)*60
              RA='00'[len(str(ra0)):]+str(ra0)+':'+'00'[len(str(ra1)):]+str(ra1)+':'+'00'[len(str(int(ra2))):]+str(ra2)
      if ra and dec:          return RA, DEC
      else:                   return RA or DEC

###########################################################################

def lscastroloop(imglist,catalogue,_interactive,number1,number2,number3,_fitgeo,_tollerance1,_tollerance2,sexvec='',_guess=False,_numin=4,method='vizir',xshift=0,yshift=0):
    import floacqastrodef
    from floacqastrodef import lscastrometry2
    from numpy import median, array
    import math
    import datetime
    import time
    _imex=False
    for img in imglist:
        hdr=floacqastrodef.readhdr(img)
        _instrume=floacqastrodef.readkey3(hdr,'instrume')
    if catalogue=='inst' and _instrume in ['FLOYDSAG-kb42','FLOYDSAG-kb37','37']:
        angle=readkey3(hdr,'ROLLERDR')#posang)
        theta=(angle*pi/180.)
        CDELT0=0.000129722   # 1.3042840792028E-4   8.43604528922325E-5  #6.6888889999999995e-05
        CD1_1=(-1)*CDELT0*cos(theta)
        CD2_2=CDELT0*cos(theta)
        CD1_2=abs(CDELT0)*(abs(CDELT0)/CDELT0)*sin(theta)
        CD2_1=abs(CDELT0)*(abs(CDELT0)/CDELT0)*sin(theta)
        if not CRPIX1:        CRPIX1= readkey3(hdr,'ROTCENTX')
        else: CRPIX1= 1000.+CRPIX1
        if not CRPIX2:        CRPIX2= readkey3(hdr,'ROTCENTY')
        else: CRPIX2= 1000.+CRPIX2
        CDELT1=2
        CDELT2=2
    elif _instrume in ['fs01','fs02']:
        angle=readkey3(hdr,'ROTSKYPA')#posang)
        theta=(angle*pi/180.)
 #       pixscale=0.30*_CCDXBIN
 #       CDELT0=pixscale/3600.
        CDELT0=0.000083568667  # 8.43604528922325E-5  #6.6888889999999995e-05
        CD1_1=(-1)*CDELT0*cos(theta)
        CD2_2=CDELT0*cos(theta)
        CD1_2=abs(CDELT0)*(abs(CDELT0)/CDELT0)*sin(theta)
        CD2_1=abs(CDELT0)*(abs(CDELT0)/CDELT0)*sin(theta)
        if not CRPIX1:        CRPIX1= 1024.
        else: CRPIX1= 1024.+CRPIX1
        if not CRPIX2:        CRPIX2= 1024.
        else: CRPIX2= 1024.+CRPIX2
        CDELT1=0.000083568694*(-1)
        CDELT2=0.000083568694
    elif _instrume in ['fs03']:
        angle=readkey3(hdr,'ROTSKYPA')#posang)
        theta=(angle*pi/180.)
 #       pixscale=0.30*_CCDXBIN
 #       CDELT0=pixscale/3600.
        CDELT0=0.000083705976   #8.43604528922325E-5  #6.6888889999999995e-05
        CD1_1=(-1)*CDELT0*cos(theta)
        CD2_2=CDELT0*cos(theta)
        CD1_2=abs(CDELT0)*(abs(CDELT0)/CDELT0)*sin(theta)
        CD2_1=abs(CDELT0)*(abs(CDELT0)/CDELT0)*sin(theta)
        if not CRPIX1:        CRPIX1= 1024.
        else: CRPIX1= 1024.+CRPIX1
        if not CRPIX2:        CRPIX2= 1024.
        else: CRPIX2= 1024.+CRPIX2
        CDELT1=0.000083705976*(-1)
        CDELT2=0.000083705976
    elif _instrume.lower() in ['em03','em01']:
        theta=(angle*pi/180.)
        CDELT0=7.63077724258886e-05 #7.7361111111111123e-05 
        CD1_1=(-1)*CDELT0*cos(theta)
        CD2_2=CDELT0*cos(theta)
        CD1_2=abs(CDELT0)*(abs(CDELT0)/CDELT0)*sin(theta)
        CD2_1=abs(CDELT0)*(abs(CDELT0)/CDELT0)*sin(theta)
        if not CRPIX1:        CRPIX1= readkey3(hdr,'ROTCENTX')
        if not CRPIX2:        CRPIX2= readkey3(hdr,'ROTCENTY')-100
        CDELT1=2
        CDELT2=2
    else:  print '\n### ERROR: instument not found !!!'

    CTYPE1  = 'RA---TAN'  
    CTYPE2  = 'DEC--TAN' 
    CRVAL1=_RA
    CRVAL2=_DEC
    WCSDIM  =                   2  
    LTM1_1  =                   1. 
    LTM2_2  =                   1.
    WAT0_001= 'system=image'
    WAT1_001= 'wtype=tan axtype=ra'
    WAT2_001= 'wtype=tan axtype=dec'
    updateheader(img,0,{'CTYPE1':[CTYPE1,''], 'CTYPE2':[CTYPE2,''], 'CRVAL1':[CRVAL1,''], 'CRVAL2':[CRVAL2,''],\
                'CRPIX1':[CRPIX1,''], 'CRPIX2':[CRPIX2,''], 'CDELT1':[CDELT1,''], 'CDELT2':[CDELT2,''],\
                'CD1_1':[CD1_1,''], 'CD2_2':[CD2_2,''], 'CD1_2':[CD1_2,''], 'CD2_1':[CD2_1,''],\
                'WCSDIM':[WCSDIM,'']})

##########################################################

def sextractor(img):
        import os
        from pyraf import iraf
        from iraf import proto
        from numpy import compress,array,asarray
        
        namesex ='default1.sex'
        os.system('sex '+img+' -c '+namesex+' > _logsex')
        delete('_logsex')
        xpix=iraf.proto.fields('detections.cat',fields='2',Stdout=1)
        ypix=iraf.proto.fields('detections.cat',fields='3',Stdout=1)
        cm=iraf.proto.fields('detections.cat',fields='4',Stdout=1)
        cl=iraf.proto.fields('detections.cat',fields='7',Stdout=1)
        fw=iraf.proto.fields('detections.cat',fields='8',Stdout=1)
        ell=iraf.proto.fields('detections.cat',fields='9',Stdout=1)
        bkg=iraf.proto.fields('detections.cat',fields='10',Stdout=1)

        cl=compress((array(xpix)!=''),array(cl,float))
        cm=compress((array(xpix)!=''),array(cm,float))
        fw=compress((array(xpix)!=''),array(fw,float))
        ell=compress((array(xpix)!=''),array(ell,float))
        bkg=compress((array(xpix)!=''),array(bkg,float))
        ypix=compress((array(xpix)!=''),array(ypix,float))
        xpix=compress((array(xpix)!=''),array(xpix,float))
        try:
            ww=asarray([i for i in range(len(xpix)) if ((xpix[i]<2010) or (ypix[i]<2010))])
            cl,cm,fw,ell,xpix,ypix,bkg=cl[ww],cm[ww],fw[ww],ell[ww],xpix[ww],ypix[ww],bkg[ww]

            ww=asarray([i for i in range(len(xpix)) if ((xpix[i]>20) or (ypix[i]<2010))])
            cl,cm,fw,ell,xpix,ypix,bkg=cl[ww],cm[ww],fw[ww],ell[ww],xpix[ww],ypix[ww],bkg[ww]

            ww=asarray([i for i in range(len(xpix)) if (xpix[i]>3)])
            cl,cm,fw,ell,xpix,ypix,bkg=cl[ww],cm[ww],fw[ww],ell[ww],xpix[ww],ypix[ww],bkg[ww]

            cl=compress((array(fw)<=15)&(array(fw)>=-2),array(cl))
            cm=compress((array(fw)<=15)&(array(fw)>=-2),array(cm))
            xpix=compress((array(fw)<=15)&(array(fw)>=-2),array(xpix))
            ypix=compress((array(fw)<=15)&(array(fw)>=-2),array(ypix))
            ell=compress((array(fw)<=15)&(array(fw)>=-2),array(ell))
            bkg=compress((array(fw)<=15)&(array(fw)>=-2),array(bkg))
            fw=compress((array(fw)<=15)&(array(fw)>=-2),array(fw))
        except: 
            xpix,ypix,fw,cl,cm,ell=[],[],[],[],[],[]
            print '\n### ERROR Filtering the sextractor detections, please check that sextractor is working ......'
        delete('detections.cat')
        return xpix,ypix,fw,cl,cm,ell,bkg

###############################################################################

def defsex(namefile):
    import lsc
    import string,re,os
    sexfile=lsc.__path__[0]+'/standard/sex/default.sex'
    f=open(sexfile,'r')
    ss=f.readlines()
    f.close()   
    ff=open(namefile,'w')
    for i in ss:
        if string.count(i,'PARAMETERS_NAME')==1:
            ff.write('PARAMETERS_NAME  "'+lsc.__path__[0]+'/standard/sex/default.param"\n')
        elif string.count(i,'FILTER_NAME')==1:
            ff.write('FILTER_NAME  "'+lsc.__path__[0]+'/standard/sex/default.conv"\n')
        elif string.count(i,'STARNNW_NAME')==1:
            ff.write('STARNNW_NAME "'+lsc.__path__[0]+'/standard/sex/default.nnw"\n')
        else:
            ff.write(i)
    ff.close()
    return namefile

###############################################################################################

def vizq(_ra,_dec,catalogue,radius):
    ''' Query vizquery '''
    import os,string,re
    import floacqastrodef
    from numpy import array,compress
#    _site='vizier.cfa.harvard.edu'
    _site='vizier.u-strasbg.fr'
    cat={'usnoa2':['I/252/out','USNO-A2.0','USNO-A2.0,Rmag'],\
         '2mass':['II/246/out','2MASS','2MASS,Jmag'],\
         'landolt':['II/183A/table2','','Vmag,B-V,U-B,V-R,R-I,Star,e_Vmag'],\
         'apass':['I/322A/out','','Bmag,Vmag,gmag,rmag,imag,e_Vmag,e_Bmag,e_gmag,e_rmag,e_imag,UCAC4'],\
         'usnob1':['I/284/out','USNO-B1.0','R2mag'],'sdss7':['II/294/sdss7','','objID,umag,gmag,rmag,imag,zmag,gc'],\
         'sdss9':['V/139/sdss9','','objID,umag,gmag,rmag,imag,zmag,e_umag,e_gmag,e_rmag,e_imag,e_zmag,gc'],\
         'sdss7':['II/294/sdss7','','objID,umag,gmag,rmag,imag,zmag,e_umag,e_gmag,e_rmag,e_imag,e_zmag,gc'],\
         'sdss8':['II/306/sdss8','','objID,umag,gmag,rmag,imag,zmag,e_umag,e_gmag,e_rmag,e_imag,e_zmag,gc']}

    a=os.popen('vizquery -mime=tsv  -site='+_site+' -source='+cat[catalogue][0]+\
                   ' -c.ra='+str(_ra)+' -c.dec='+str(_dec)+' -c.eq=J2000 -c.rm='+str(radius)+\
                   ' -c.geom=b -oc.form=h -sort=_RA*-c.eq -out.add=_RAJ2000,_DEJ2000 -out.max=10000 -out='+\
                   cat[catalogue][1]+' -out='+cat[catalogue][2]+'').read()
    print 'vizquery -mime=tsv  -site='+_site+' -source='+cat[catalogue][0]+\
                   ' -c.ra='+str(_ra)+' -c.dec='+str(_dec)+' -c.eq=J2000 -c.rm='+str(radius)+\
                   ' -c.geom=b -oc.form=h -sort=_RA*-c.eq -out.add=_RAJ2000,_DEJ2000 -out.max=10000 -out='+\
                   cat[catalogue][1]+' -out='+cat[catalogue][2]+''
    aa=a.split('\n')
    bb=[]
    for i in aa:
        if i and i[0]!='#':   bb.append(i)
    _ra,_dec,_name,_mag=[],[],[],[]
    for ii in bb[3:]:
        aa=ii.split('\t')
        rr,dd=floacqastrodef.deg2HMS(ra=re.sub(' ',':',aa[0]), dec=re.sub(' ',':',aa[1]), round=False)
        _ra.append(rr)
        _dec.append(dd)
        _name.append(aa[2])
    dictionary={'ra':_ra,'dec':_dec,'id':_name}
    sss=string.split(cat[catalogue][2],',')
    for ii in sss: dictionary[ii]=[]
    for ii in bb[3:]:
        aa=ii.split('\t')
        for gg in range(0,len(sss)):
           if sss[gg] not in ['UCAC4','id']:
              try:
                 dictionary[sss[gg]].append(float(aa[2+gg]))
              except:    
                 dictionary[sss[gg]].append(float(9999))
           else:
                 dictionary[sss[gg]].append(str(aa[2+gg]))

    if catalogue in ['sdss7','sdss9','sdss8']:
        dictionary['u']=dictionary['umag']
        dictionary['g']=dictionary['gmag']
        dictionary['r']=dictionary['rmag']
        dictionary['i']=dictionary['imag']
        dictionary['z']=dictionary['zmag']
        dictionary['uerr']=dictionary['e_umag']
        dictionary['gerr']=dictionary['e_gmag']
        dictionary['rerr']=dictionary['e_rmag']
        dictionary['ierr']=dictionary['e_imag']
        dictionary['zerr']=dictionary['e_zmag']
        for key in dictionary.keys():
           if key!='r':
              dictionary[key]=compress((array(dictionary['r'])<19)&(array(dictionary['r']>13)),dictionary[key])
        dictionary['r']=compress((array(dictionary['r'])<19)&(array(dictionary['r']>13)),dictionary['r'])

    elif  catalogue=='landolt':
        dictionary['B']=array(dictionary['Vmag'])+array(dictionary['B-V'])
        dictionary['U']=array(dictionary['B'])+array(dictionary['U-B'])
        dictionary['V']=array(dictionary['Vmag'])
        dictionary['Verr']=array(dictionary['e_Vmag'])
        dictionary['R']=array(dictionary['Vmag'])-array(dictionary['V-R'])
        dictionary['I']=array(dictionary['R'])-array(dictionary['R-I'])
        dictionary['id']=array(dictionary['Star'])
    elif  catalogue=='apass':
        dictionary['B']=array(dictionary['Bmag'])
        dictionary['V']=array(dictionary['Vmag'])
        dictionary['g']=array(dictionary['gmag'])
        dictionary['r']=array(dictionary['rmag'])
        dictionary['i']=array(dictionary['imag'])
        dictionary['Berr']=array(dictionary['e_Bmag'],float)/100.
        dictionary['Verr']=array(dictionary['e_Vmag'],float)/100.
        dictionary['gerr']=array(dictionary['e_gmag'],float)/100.
        dictionary['rerr']=array(dictionary['e_rmag'],float)/100.
        dictionary['ierr']=array(dictionary['e_imag'],float)/100.
        dictionary['id']=array(dictionary['UCAC4'],str)
        for key in dictionary.keys():
           if key!='r':
              dictionary[key]=compress((array(dictionary['r'])<19)&(array(dictionary['r']>10.5)),dictionary[key])
        dictionary['r']=compress((array(dictionary['r'])<19)&(array(dictionary['r']>10.5)),dictionary['r'])
    return dictionary
####################################################################################################

##############################################################################
def delete(listfile):
    import os,string,re,glob
    if listfile[0]=='@':   
        ff = open(listfile[1:])
        files = ff.readlines()
        imglist = []
        for ff in files: 
            ff=re.sub(' ','',ff)
            if not ff=='\n' and ff[0]!='#':
                ff=re.sub('\n','',ff)
                imglist.append(ff)
    elif ',' in listfile: imglist = string.split(listfile,sep=',')
    else:       imglist=[listfile]    
    lista=[]
    for _file in imglist:   lista=lista+glob.glob(_file)
    if lista:
        for _file in lista:
            try:          os.system('rm '+_file)
            except:       pass
###############################################################
#############################################

def deg2HMS(ra='', dec='', round=False):
      import string
      RA, DEC= '', ''
      if dec:
          if string.count(str(dec),':')==2:
              dec00=string.split(dec,':')
              dec0,dec1,dec2=float(dec00[0]),float(dec00[1]),float(dec00[2])
              if '-' in str(dec0):       DEC=(-1)*((dec2/60.+dec1)/60.+((-1)*dec0))
              else:                      DEC=(dec2/60.+dec1)/60.+dec0
          else:
              if str(dec)[0]=='-':      dec0=(-1)*abs(int(dec))
              else:                     dec0=abs(int(dec))
              dec1=int((abs(dec)-abs(dec0))*(60))
              dec2=((((abs(dec))-abs(dec0))*60)-abs(dec1))*60
              DEC='00'[len(str(dec0)):]+str(dec0)+':'+'00'[len(str(dec1)):]+str(dec1)+':'+'00'[len(str(int(dec2))):]+str(dec2)
      if ra:
          if string.count(str(ra),':')==2:
              ra00=string.split(ra,':')
              ra0,ra1,ra2=float(ra00[0]),float(ra00[1]),float(ra00[2])
              RA=((ra2/60.+ra1)/60.+ra0)*15.
          else:
              ra0=int(ra/15.)
              ra1=int(((ra/15.)-ra0)*(60))
              ra2=((((ra/15.)-ra0)*60)-ra1)*60
              RA='00'[len(str(ra0)):]+str(ra0)+':'+'00'[len(str(ra1)):]+str(ra1)+':'+'00'[len(str(int(ra2))):]+str(ra2)
      if ra and dec:          return RA, DEC
      else:                   return RA or DEC

###########################################################################

def lscastroloop(imglist,catalogue,_interactive,number1,number2,number3,_fitgeo,_tollerance1,_tollerance2,sexvec='',_guess=False,_numin=4,method='vizir',xshift=0,yshift=0):
    import floacqastrodef
    from floacqastrodef import lscastrometry2
    from numpy import median, array
    import math
    import datetime
    import time
    _imex=False
    for img in imglist:
        hdr=floacqastrodef.readhdr(img)
        _instrume=floacqastrodef.readkey3(hdr,'instrume')
        if catalogue=='inst' and _instrume in ['FLOYDSAG-kb42','FLOYDSAG-kb37','kb37']: catalogue='2mass'
        if not sexvec:
            sexvec=floacqastrodef.sextractor(img)
###################
        print xshift,yshift
        if xshift!=0 and yshift!=0:
            print 'guess astrometry before starting '
            floacqastrodef.wcsstart(img,xshift,yshift)
        catvec=floacqastrodef.querycatalogue(catalogue,img,method)
        #print catalogue,sexvec,catvec
        rmsx1,rmsy1,num1,fwhm1,ell1,ccc,bkg1,rasys1,decsys1=floacqastrodef.lscastrometry2([img],catalogue,_interactive,number1,sexvec,catvec,guess=False,fitgeo=_fitgeo,\
                                                                                 tollerance1=_tollerance1, tollerance2=_tollerance2,_update='yes',imex=_imex,nummin=_numin)
        if rmsx1>1 or rmsy1>1:
            catvec=floacqastrodef.querycatalogue(catalogue,img,method)
            rmsx2,rmsy2,num2,fwhm2,ell2,ccc,bkg2,rasys2,decsys2=floacqastrodef.lscastrometry2([img],catalogue,_interactive,number2,sexvec,catvec,guess=False,fitgeo=_fitgeo,\
                                                                                     tollerance1=_tollerance1, tollerance2=_tollerance2,_update='yes',imex=_imex,nummin=_numin)
            if rmsx2>1 or rmsy2>1:
                catvec=floacqastrodef.querycatalogue(catalogue,img,method)
                rmsx3,rmsy3,num3,fwhm3,ell3,ccc,bkg3,rasys3,decsys3=floacqastrodef.lscastrometry2([img],catalogue,_interactive,number3,sexvec,catvec,guess=False,fitgeo=_fitgeo,\
                                                                                         tollerance1=_tollerance1, tollerance2=_tollerance2,_update='yes',imex=_imex,nummin=_numin)
            else:  rmsx3,rmsy3,num3,fwhm3,ell3,ccc,bkg3,rasys3,decsys3=rmsx2,rmsy2,num2,fwhm2,ell2,ccc,bkg2,rasys2,decsys2
        else:  rmsx3,rmsy3,num3,fwhm3,ell3,ccc,bkg3,rasys3,decsys3=rmsx1,rmsy1,num1,fwhm1,ell1,ccc,bkg1,rasys1,decsys1
######################################## 
        if rmsx3 < 10 and rmsy3 < 10: 
            if _instrume in ['FLOYDSAG-kb42','FLOYDSAG-kb37','kb37']:
                fwhmgess3=median(array(fwhm3))*.68*2.35*0.467
                if _imex:  fwhmgessime = median(array(ccc))*0.467
                else:     fwhmgessime = 9999
            elif _instrume in ['fs01','fs02','fs03']:
                fwhmgess3=median(array(fwhm3))*.68*2.35*0.30
                if _imex:  fwhmgessime = median(array(ccc))*0.30
                else:     fwhmgessime = 9999
            elif _instrume in ['em03','em01']:
                fwhmgess3=median(array(fwhm3))*.68*2.35*0.278
                if _imex:  fwhmgessime = median(array(ccc))*0.278  
                else:     fwhmgessime = 9999
            ellgess3=median(array(ell3))
        else:
            fwhmgess3=9999
            fwhmgessime=9999
            ellgess3=9999
        if _instrume in ['FLOYDSAG-kb42','FLOYDSAG-kb37','kb37']:
            mbkg3=median(bkg3)
            floacqastrodef.updateheader(img,0,{'MBKG':[mbkg3,'background level']})
        else:
            mbkg3=floacqastrodef.readkey3(hdr,'MBKG')
        if fwhmgess3:
            print fwhmgess3
            if _instrume in ['FLOYDSAG-kb42','FLOYDSAG-kb37','kb37']:
                V=(math.pi/(4*math.log(2)))*(45000-float(mbkg3))*(float(fwhmgess3)**2)
            else:                     
                V=(math.pi/(4*math.log(2)))*(32000-float(mbkg3))*(float(fwhmgess3)**2)
            magsat=-2.5*math.log10(V)
        else:        magsat=9999
    print rmsx3,rmsy3,num3,fwhmgess3,ellgess3,fwhmgessime,rasys3,decsys3,magsat
    return  rmsx3,rmsy3,num3,fwhmgess3,ellgess3,fwhmgessime,rasys3,decsys3,magsat

##################################################################################################
############################################################################################################
def lscastrometry2(lista,catalogue,_interactive,number,sexvec,catvec,guess=False,fitgeo='xyscale', tollerance1=100, tollerance2=30, _update='yes',imex=False,nummin=4):
    import os,string,re,sys
    import numpy
    import math
    from numpy import array, compress, argsort, sort, asarray
    from numpy import round, mean, std, sqrt, median
    from numpy import argmin, isnan, abs, genfromtxt
    import time
    import datetime
    import floacqastrodef
    from floacqastrodef import wcsstart
    from floacqastrodef import delete, readhdr,readkey3, display_image
    from pyraf import iraf
    xpix,ypix,fw,cl,cm,ell,bkg=sexvec
    acoo1,apix1,am1=catvec['coo'],catvec['pix'],catvec['mag']
########################   catalogue
    iraf.noao(_doprint=0)
    iraf.imcoords(_doprint=0)
    iraf.tv(_doprint=0)
    iraf.tv.rimexam.backgrou = 'yes'
    iraf.astcat(_doprint=0)
    toforget = ['imcoords','astcat','tv']
    for t in toforget: iraf.unlearn(t)
    verbose=False
    if _interactive: verbose=True
    img=lista[0]
    hdr=readhdr(img)
    _instrume=readkey3(hdr,'instrume')
    if _instrume in ['FLOYDSAG-kb42','FLOYDSAG-kb37','kb37']:
            magsel0 = 7.0
            magsel1 = 21.
    else:
            magsel0 = 7.0
            magsel1 = 21.
    _CRPIX1=readkey3(hdr,'CRPIX1')
    _CRPIX2=readkey3(hdr,'CRPIX2')
    if verbose:            display_image(img,1,'','',False)
    if verbose:
            iraf.tvmark(1,'STDIN',Stdin=list(apix1),mark="circle",number='no',label='no',radii=5,nxoffse=5,nyoffse=5,color=205,txsize=2)
            raw_input('mark catalogue '+str(len(apix1)))
    else:  
#        ss=datetime.datetime.now()
        time.sleep(.7)
    answ='yes'
    magsel11=magsel1
    mlim=0
    while answ=='yes':
            amcut1=compress((array(am1)>magsel0) &(array(am1)<magsel11), am1)
            if len(amcut1)<=number:    
                answ='no'
                magsel11=magsel1+mlim+.5   
            else:
                mlim=mlim-.5
                magsel11=magsel1+mlim

    amcut=compress((array(am1)>magsel0) &(array(am1)<magsel11), am1)
    apixcut=compress((array(am1)>magsel0) &(array(am1)<magsel11), apix1)      #   usno x y  cut_list
    acoocut=compress((array(am1)>magsel0) &(array(am1)<magsel11), acoo1)  #   usno ra dec  cut_list

    rausno=compress((array(am1)>magsel0) &(array(am1)<magsel11), array(catvec['ra'],float)) 
    decusno=compress((array(am1)>magsel0) &(array(am1)<magsel11), array(catvec['dec'],float)) 
    xusno,yusno=[],[]
    for i in apixcut:
            xusno.append(float(string.split(i)[0]))
            yusno.append(float(string.split(i)[1]))
    xusno,yusno=array(xusno),array(yusno)
    
#################################################################
    if verbose:
            iraf.tvmark(1,'STDIN',Stdin=list(apixcut),mark="circle",number='yes',label='no',radii=8,nxoffse=5,nyoffse=5,color=204,txsize=2)
            raw_input('brightest '+str(number)+' objects')

##############    sextractor   ##################
    if len(xpix)>=number:
            cm=array(cm,float)
            xpix=xpix[argsort(cm)][0:number]
            ypix=ypix[argsort(cm)][0:number]
            fw=fw[argsort(cm)][0:number]
            ell=ell[argsort(cm)][0:number]
            cm=cm[argsort(cm)][0:number]
    if verbose:
            sexpix=[]
            for i in range(0,len(xpix)):
                sexpix.append(str(xpix[i])+' '+str(ypix[i]))
            iraf.tvmark(1,'STDIN',Stdin=list(sexpix),mark="circle",number='no',label='no',radii=8,nxoffse=5,nyoffse=5,color=206,txsize=2)
            raw_input('print sex '+str(len(sexpix)))

    xsex,ysex=array(xpix),array(ypix)
    fwsex=array(fw)
    ellsex=array(ell)
#####################################################################
    max_sep=tollerance1
    xdist,ydist=[],[]
    for i in range(len(xusno)):
            dist = sqrt((xusno[i]-xsex)**2+(yusno[i]-ysex)**2)
            idist = argmin(dist)
            if dist[idist]<max_sep:
                xdist.append(xusno[i]-xsex[idist])
                ydist.append(yusno[i]-ysex[idist])
    xoff,xstd = round(median(xdist),2),round(std(xdist),2)
    yoff,ystd = round(median(ydist),2),round(std(ydist),2)
    _xdist,_ydist = array(xdist),array(ydist)
    __xdist = compress((abs(_xdist-xoff)<3*xstd)&(abs(_ydist-yoff)<3*ystd),_xdist)
    __ydist = compress((abs(_xdist-xoff)<3*xstd)&(abs(_ydist-yoff)<3*ystd),_ydist)
    xoff,xstd = round(median(__xdist),2),round(std(__xdist),2)
    yoff,ystd = round(median(__ydist),2),round(std(__ydist),2)
    if isnan(xoff): xoff=0
    if isnan(yoff): yoff=0
    _CRPIX1=readkey3(hdr,'CRPIX1')
    _CRPIX2=readkey3(hdr,'CRPIX2')
    floacqastrodef.updateheader(img,0,{'CRPIX1':[_CRPIX1-xoff,'']})
    floacqastrodef.updateheader(img,0,{'CRPIX2':[_CRPIX2-yoff,'']})
    xusno2_new=xusno-xoff
    yusno2_new=yusno-yoff
#####################################################################
    max_sep=tollerance2
    fwhm=[]
    fwhm2=[]
    ell=[]
    xref=[]
    iraf.tv(_doprint=0)
    iraf.tv.rimexam.backgrou = 'yes'
    vettoretran=[]
    for i in range(len(xusno2_new)):
            dist = sqrt((xusno2_new[i]-xsex)**2+(yusno2_new[i]-ysex)**2)
            idist = argmin(dist)
            if dist[idist]<max_sep:
                xref.append(xsex[idist])
                vettoretran.append(str(rausno[i])+' '+str(decusno[i])+' '+str(xsex[idist])+' '+str(ysex[idist])+' \n')
                fwhm.append(fwsex[idist])
                ell.append(ellsex[idist])
                if imex:
                    gg=open('tmp.one','w')
                    gg.write(str(xsex[idist])+' '+str(ysex[idist])+'\n')
                    gg.close()
                    ime=iraf.imexam(input=img, frame=1, logfile='', keeplog='yes', imagecur='tmp.one', wcs='logical', use_disp='no',Stdout=1)
                    try:
                        _fwhm2=median(compress(array(string.split(ime[3])[-3:],float)<99,(array(string.split(ime[3])[-3:],float))))
                        fwhm2.append(_fwhm2)
                    except: pass
    if len(xref)>=nummin:
            _ccmap1=iraf.ccmap('STDIN','STDOUT',images=img,Stdin=vettoretran,fitgeome=fitgeo,xcolum=3, xxorder=2,\
                               yyorder=2, ycolum=4,lngcolum=1,latcolumn=2,lngunit='degrees',update='No',interact='No',maxiter=3,Stdout=1)
            if 'rms' in _ccmap1[_ccmap1.index('Wcs mapping status')+1]:
                try:           rmsx,rmsy=array(string.split(string.split(_ccmap1[_ccmap1.index('Wcs mapping status')+1],':')[-1])[0:2],float)
                except:        rmsx,rmsy=array(string.split(string.split(_ccmap1[_ccmap1.index('Wcs mapping status')+1],':')[-1])[0:2])
                if rmsx<2 and rmsy<2:
                    print '\n### update astrometry with non linear order' 
                    _ccmap1=iraf.ccmap('STDIN','STDOUT',images=img,Stdin=vettoretran,fitgeome=fitgeo,xcolum=3, xxorder=2,\
                                       yyorder=2, ycolum=4,lngcolum=1,latcolumn=2,lngunit='degrees',update='Yes',interact='No',maxiter=3,Stdout=1)
                    xy = iraf.wcsctran('STDIN',output="STDOUT",Stdin=vettoretran,Stdout=1,image=img,inwcs='physical', outwcs='world',column="3 4",formats='%10.6f %10.6f',verbose='yes')[3:]
                    rasys=median(array([float(xy[i].split()[2])-float(xy[i].split()[0]) for i in range(0, len(xy)) ]))
                    decsys=median(array([float(xy[i].split()[3])-float(xy[i].split()[1]) for i in range(0, len(xy)) ]))
                else:
                    rasys,decsys=999,999
            else:      rmsx,rmsy,rasys,decsys=999,999,999,999
    else:
            rmsx,rmsy=990,999
            rasys,decsys=999,999
    return rmsx,rmsy,len(xref),fwhm,ell,fwhm2,bkg,rasys,decsys
##########################################################################


def querycatalogue(catalogue,img,method='vizir'):
        from pyraf import iraf
        from numpy import array, compress
        import string,re,sys,os
        import floacqastrodef
        from floacqastrodef import delete,readhdr,readkey3
        hdr=readhdr(img)
        _ra=readkey3(hdr,'RA')
        _dec=readkey3(hdr,'DEC')
        _ra,_dec=floacqastrodef.deg2HMS(_ra,_dec)
        stdcoo={}
        if method=='vizir':
            stdcoo=floacqastrodef.vizq(_ra,_dec,catalogue,20)
            lll=['# END CATALOG HEADER','#']
            colonne4={'usnoa2':'Rmag','usnob1':'mag','2mass':'Jmag','gsc1':'mag'}
            for ff in range(0,len(stdcoo['ra'])):
                lll.append(str(stdcoo['ra'][ff])+'  '+str(stdcoo['dec'][ff])+'  '+str(stdcoo[colonne4[catalogue]][ff]))
            #print lll
            colonne3=' 1   2 '
            column={'ra':1,'dec':2,'r':3}
        if  string.count(str(stdcoo['ra'][0]),':'):
            ddd2=iraf.wcsctran('STDIN','STDOUT',img,Stdin=lll,Stdout=1,inwcs='world',units='hour degrees',outwcs='logical',columns=colonne3,formats='%10.1f %10.1f')
        else:
            ddd2=iraf.wcsctran('STDIN','STDOUT',img,Stdin=lll,Stdout=1,inwcs='world',units='degree degrees',outwcs='logical',columns=colonne3,formats='%10.1f %10.1f')
        xx,yy=[],[]
        for i in ddd2[ddd2.index('# END CATALOG HEADER')+2:]:
            xx.append(float(i.split()[column['ra']-1]))
            yy.append(float(i.split()[column['dec']-1]))
#######
        acoo1=[]
        apixx1,apixy1,am1,apix1=[],[],[],[]
        for i in range(0,len(stdcoo['ra'])):
            acoo1.append(str(stdcoo['ra'][i])+' '+str(stdcoo['dec'][i]))
            apix1.append(str(xx[i])+' '+str(yy[i]))
            am1.append(stdcoo[colonne4[catalogue]][i])
            if  string.count(str(stdcoo['ra'][i]),':'):
                stdcoo['ra'][i]=(int(string.split(stdcoo['ra'][i],':')[0])+float(string.split(stdcoo['ra'][i],':')[1])/60+float(string.split(stdcoo['ra'][i],':')[2])/3600.)*15
                if string.count(str(stdcoo['dec'][i]),'-')==0:   
                    stdcoo['dec'][i]=int(string.split(stdcoo['dec'][i],':')[0])+float(string.split(stdcoo['dec'][i],':')[1])/60+float(string.split(stdcoo['dec'][i],':')[2])/3600.
                else:
                    stdcoo['dec'][i]=(-1)*(abs(int(string.split(stdcoo['dec'][i],':')[0]))+float(string.split(stdcoo['dec'][i],':')[1])/60+float(string.split(stdcoo['dec'][i],':')[2])/3600.)
        stdcoo['ra']=array(stdcoo['ra'],float)
        stdcoo['dec']=array(stdcoo['dec'],float)
        if catalogue=='2mass':
            for jj in range(0,len(am1)):
                    try: am1[jj]=float(re.sub('L','',str(am1[jj])))
                    except: am1[jj]=999

        for key in stdcoo.keys():
            try:
                stdcoo[key]=compress((array(xx)<int(int(hdr['NAXIS1'])+100))&(array(xx)>-100)&(array(yy)<int(int(hdr['NAXIS2'])+100))&(array(yy)>-100),array(stdcoo[key]))
            except:  pass
        stdcoo['coo']=compress((array(xx)<int(int(hdr['NAXIS1'])+100))&(array(xx)>-100)&(array(yy)<int(int(hdr['NAXIS2'])+100))&(array(yy)>-100),array(acoo1))
        stdcoo['pix']=compress((array(xx)<int(int(hdr['NAXIS1'])+100))&(array(xx)>-100)&(array(yy)<int(int(hdr['NAXIS2'])+100))&(array(yy)>-100),array(apix1))
        stdcoo['mag']=compress((array(xx)<int(int(hdr['NAXIS1'])+100))&(array(xx)>-100)&(array(yy)<int(int(hdr['NAXIS2'])+100))&(array(yy)>-100),array(am1,float))
        stdcoo['x']=compress((array(xx)<int(int(hdr['NAXIS1'])+100))&(array(xx)>-100)&(array(yy)<int(int(hdr['NAXIS2'])+100))&(array(yy)>-100),array(xx,float))
        stdcoo['y']=compress((array(xx)<int(int(hdr['NAXIS1'])+100))&(array(xx)>-100)&(array(yy)<int(int(hdr['NAXIS2'])+100))&(array(yy)>-100),array(yy,float))
        return stdcoo

########################################################################################################################################
#################################################################################################
def display_image(img,frame,_z1,_z2,scale,_xcen=0.5,_ycen=0.5,_xsize=1,_ysize=1,_erase='yes'):
    goon='True'
    import glob, subprocess, os, time
    ds9 = subprocess.Popen("ps -U"+str(os.getuid())+"|grep -v grep | grep ds9",shell=True,stdout=subprocess.PIPE).stdout.readlines()
    if len(ds9)== 0 :   
       subproc = subprocess.Popen('ds9',shell=True)
       time.sleep(3)

    if glob.glob(img):
       from pyraf import iraf
       iraf.images(_doprint=0)
       iraf.tv(_doprint=0)
       import string,os
       if _z2: 
          try:
              sss=iraf.display(img, frame, xcen=_xcen, ycen=_ycen, xsize=_xsize, ysize=_ysize, erase=_erase,\
                                   fill='yes', zscale='no', zrange='no', z1=_z1, z2=_z2,Stdout=1)
          except:
              print ''
              print '### ERROR: PROBLEM OPENING DS9'
              print ''
              goon='False'                 
       else:
        try:  
            sss=iraf.display(img, frame, xcen=_xcen, ycen=_ycen, xsize=_xsize, ysize=_ysize, erase=_erase, fill='yes', Stdout=1)
        except:
            print ''
            print '### ERROR: PROBLEM OPENING DS9'
            print ''
            goon=False
 
       if scale and goon:
          answ0 = raw_input('>>> Cuts OK ? [y/n] ? [y] ')
          if not answ0: answ0='y'
          elif answ0=='no' or answ0=='NO': answ0='n' 

          while answ0=='n':
              _z11=float(string.split(string.split(sss[0])[0],'=')[1])
              _z22=float(string.split(string.split(sss[0])[1],'=')[1])
              z11 = raw_input('>>> z1 = ? ['+str(_z11)+'] ? ')
              z22 = raw_input('>>> z2 = ? ['+str(_z22)+'] ? ')
              if not z11: z11=_z11
              else: z11=float(z11)
              if not z22: z22=_z22
              else: z22=float(z22)
              print z11,z22
              sss=iraf.display(img,frame,fill='yes', xcen=_xcen, ycen=_ycen, xsize=_xsize, ysize=_ysize, erase=_erase,\
                                   zrange='no', zscale='no', z1=z11, z2=z22, Stdout=1)
              answ0 = raw_input('>>> Cuts OK ? [y/n] ? [y] ')
              if not answ0: answ0='y'
              elif answ0=='no' or answ0=='NO': answ0='n'
       if goon:
          _z1,_z2=string.split(string.split(sss[0])[0],'=')[1],string.split(string.split(sss[0])[1],'=')[1]
    else:
        print 'Warning: image '+str(img)+' not found in the directory '
    return _z1,_z2,goon

#################################################################################################################
#################################################################################################################

def findcenter(img,hw=30,verbose=False):
    """
    this function find the pixel position of the slit image 
    on the FLOYDS acquisition images kb37, kb42  
    is cross-corelate a delta function with counts along a column passing trough the slit image
    hw = the width of the delta maximum and should be related to the slith width    
    """
    #   this is the average slope of the slit on the acquisition image
    slope=0.109 
    def gaussian(x, mu, sig):
        return np.exp(-np.power(x - mu, 2.) / 2 * np.power(sig, 2.))
    try:
        print 'here'
        data,hdr0 = pyfits.getdata(img, header=True)   #    read image
        #
        #######   we are cutting the image close to the slit position
        #######   parameters are different depending if the image is binned or not 
        #
        if len(data)<600:
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
        data1=data[y0:y1,xzero-50*rr:xzero+50*rr]      #  cut array close to slit position 
        if verbose:
            import pylab as plt
            plt.clf()
            plt.ion()
            ax1=plt.axes([.1,.1,.4,.8])
            ax2=plt.axes([.5,.1,.4,.8])
            ax1.imshow(data1, interpolation='nearest', cmap=plt.cm.ocean)

        prewvalue=''
        xxvec=[]
        yyvec=[]
        ##################################
        for nn in range(-25,25,5):
            pixelpos=(50+nn)*rr
            xx=np.arange(len(data1[:,pixelpos]))
            xx1=np.arange(xx[0],xx[-1],.1)
            #     define column  
            yy=(data1[:,pixelpos:pixelpos+4].mean(1)-max(data1[:,pixelpos:pixelpos+4].mean(1)))*(-1)  
            yy1=np.interp(xx1,xx,yy)
            zz1=np.zeros(len(yy1))+0.001
            zz1[len(zz1)/2-hw:len(zz1)/2+hw]=np.max(yy1)
            if prewvalue:
                ff=gaussian(np.arange(len(yy1)),prewvalue, 0.01)
                xcorr=np.correlate(yy1*ff,zz1,mode='same')
            else:
                xcorr=np.correlate(yy1,zz1,mode='same')
            #  exclude the edges  
            xcorr[0:len(xcorr)/10]=xcorr.min()
            xcorr[-len(xcorr)/10:]=xcorr.min()
            if np.median(yy)-((np.max(yy)-np.min(yy))*2/3)<0:
                yyvec.append(xcorr.argmax()/10.+y0)
                xxvec.append(pixelpos+4+xzero-50*rr)
                prewvalue=xcorr.argmax()+slope*5
                if verbose:
                    ax2.cla()
                    ax2.plot(np.arange(len(xcorr)),xcorr,'r-')
#                    ax2.plot(xx1,yy1,'r-')
#                    ax2.plot(xx1,zz1,'b-')
                    plt.draw()
                    ax1.plot([pixelpos+2.],[xcorr.argmax()/10.],'ro')
                    plt.draw()
            else:
                if prewvalue:  prewvalue=prewvalue+slope*5
        #################################
        #  linear fit on the slit position
        if len(xxvec)>=4:
            aa,bb=np.polyfit(xxvec,yyvec,1)
            if max(yyvec)-min(yyvec)<15 and np.mean(yyvec)>y0+10 and np.mean(yyvec)<y1-10:
                ypos=bb+xzero*aa
                xpos=xzero
            else:
                ypos=''
                xpos=''
            print xpos,ypos
        else:
            ypos=''
            xpos=''
    except  Exception as e:
        print e
        aa,bb='',''
        ypos=''
        xpos=''
    return  xpos,ypos,aa,bb

#######################################
