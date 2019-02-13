#!/usr/bin/env python


import sys
from scipy import fftpack
from scipy import misc
import numpy as np
import struct
from numpy import *
import matplotlib.pyplot as plt

#########################################################################
# get rid of annoying future warning - need to fix this for later versions
import warnings
warnings.filterwarnings("ignore")
#########################################################################

#------------- READ MRCs--------------
# from David Stokes
class mrc_image:
    def __init__(self,filename):
        self.numbytes1=56           # 56 long ints
        self.numbytes2=80*10          # 10 header lines of 80 chars each
        self.filename=filename

def read(self):
    input_image=open(self.filename,'rb')
    self.header1=input_image.read(self.numbytes1*4)
    self.header2=input_image.read(self.numbytes2)
    byte_pattern='=' + 'l'*self.numbytes1   #'=' required to get machine independent standard size
    self.dim=struct.unpack(byte_pattern,self.header1)[:3]   #(dimx,dimy,dimz)
    self.imagetype=struct.unpack(byte_pattern,self.header1)[3]  #0: 8-bit signed, 1:16-bit signed, 2: 32-bit float, 6: unsigned 16-bit (non-std)
    if (self.imagetype == 0):
        imtype='b'
    elif (self.imagetype == 1):
        imtype='h'
    elif (self.imagetype == 2):
        imtype='f4'
    elif (self.imagetype == 6):
        imtype='H'
    else:
        imtype='unknown'   #should put a fail here
    input_image_dimension=(self.dim[1],self.dim[0])  #2D images assumed
    self.image_data=fromfile(file=input_image,dtype=imtype,count=self.dim[0]*self.dim[1]).reshape(input_image_dimension)
    input_image.close()
    return(self.image_data)
#---------------------------

#----------- calculate Power spectrum---------------------------
def get_PS(image_in,outname):
    print('\ncalculating FFT and PS')
    mrcim = mrc_image(image_in)
    image = read(mrcim)

    # fourier transform of the image.

    F1 = fftpack.fft2(image)
     
    # make pretty - put low res in center
    F2 = fftpack.fftshift( F1 )
     
    # Calculate 2D power spectrum
    ps2D = np.array(np.abs( F2 )**2,dtype=np.uint8)
    
    #save the powerspectrum
    misc.imsave('{0}_PS.tif'.format(imgname),ps2D)

    
    return (F2)
#----------------------------------------------------------------

#--- calc distance on PS---------------------------------------
def calc_dist(ps2d,pxsize,resolution):
    """ given a resolution calculate its point
    on the PS
    """
    nyquist = 2*pxsize
    distancex = ((np.shape(ps2d)[1]*.5)*nyquist)/resolution
    distancey = ((np.shape(ps2d)[0]*.5)*nyquist)/resolution

    return(distancex,distancey)
#-----------------------------------------------------------------

pixelsize = float(sys.argv[2])
nyquist = pixelsize*2.0

#make the FFT
imgname = sys.argv[1].split('/')[-1].split('.')[0]
image_FFT = get_PS(sys.argv[1],'{0}_PS.png'.format(imgname))

#find out where to cut...
cut_aty= (70,45)
cut_atx=(45,nyquist)

ycent = (0.5*image_FFT.shape[0])
xcent = (0.5*image_FFT.shape[1])

xmax = int(calc_dist(image_FFT,pixelsize,45)[0])

ymax = int(calc_dist(image_FFT,pixelsize,45)[1])
ymin = int(calc_dist(image_FFT,pixelsize,70)[1])

print('\nRetaining {0} - {1} Angstrom ({2} - {3} px) perpendicular to fibril axis'.format(cut_aty[0],cut_aty[1],ymin,ymax))
print('Discarding {0} - {1} Angstrom ({2} px to edge) parallel to fibril axis'.format(cut_atx[0],cut_atx[1],xmax))

#cut in the y direction
image_FFT[:int(ycent-ymax),0:] = 0
image_FFT[int(ycent+ymax):,0:] = 0
image_FFT[int(ycent-ymin):int(ycent+ymin),0:] = 0

# cut in the xdirection
image_FFT[0:,0:int(xcent-xmax)] = 0
image_FFT[0:,int(xcent+xmax):] = 0
inverse = fftpack.ifft2(image_FFT)

#set the pixel value range
out = np.array(np.abs(inverse),dtype=np.uint8)
factor = 255/np.max(out)
scalearray= np.full(out.shape,factor)
imarray = out*scalearray    

#save the image
misc.imsave('{0}_inv.tif'.format(imgname),imarray)

#flatten the image perpendicular to fibrilaxis
flat = np.sum(out,axis=1)
x= [i*pixelsize for i in range(len(flat))]
plt.plot(x,flat)

def moving_average(a, n) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

avgbox = 100
avg = moving_average(flat,avgbox)
avgx = [(0.5*avgbox)+(i*pixelsize) for i in range(len(avg))]
plt.plot(avgx,avg)

peaksx,peaksy= [],[]
indpeaks = []
n=0
ipcounter = -1
inpeak = False
for i in avg:
    if i < np.mean(avg)-(np.std(avg)):
        if inpeak == False:
            ipcounter +=1
            inpeak = True
        peaksy.append(i)
        peaksx.append(avgx[n])
        try:
            indpeaks[ipcounter].append(avgx[n])
        except:
            if len(indpeaks) == 0:
                indpeaks = [[avgx[n]]]
            else:
                indpeaks.append([avgx[n]])
    else:
        inpeak = False
    n+=1
plt.scatter(peaksx,peaksy,color='red')

n=1
print('\nPeaks found at (discard 1st and last): ')
for i in indpeaks[1:-1]:
    print round(np.mean(i),2)
    plt.axvline(x=np.mean(i),color='k', linestyle='--')

print('\n{0} crossovers:'.format(len(indpeaks[1:-2])))    
for i in indpeaks[1:-2]:
    print round(np.mean(indpeaks[n+1])-np.mean(i),2)
    n+=1

plt.ylabel('pixel intensity')
plt.xlabel('Length (Angstroms)')
plt.savefig('plot.png')

# save the chopped PS 
#ps2D = np.array(np.abs( image_FFT )**2,dtype=np.uint8)
#scipy.misc.imsave('inv_PS.tif'.format(imgname),ps2D)

