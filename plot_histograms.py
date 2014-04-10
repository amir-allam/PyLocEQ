import random
import matplotlib.pyplot as plt
import numpy as np
from misc_tools import *
from mtools import *

fnam='Apr7_reloc'
bins=50
a=np.loadtxt(fnam)
newlon=a[:,2]
newlat=a[:,3]
newdep=a[:,4]
oldlon=a[:,5]
oldlat=a[:,6]
olddep=a[:,7]
calctime=a[:,9]

depchange=olddep-newdep
latchange=oldlat-newlat
lonchange=oldlon-newlon
horiz_change=sqrt( (lonchange*111)**2+(latchange*111)**2 )

plt.figure("Depth");plt.clf()
plt.hist(depchange/3,bins=bins)
plt.xlim([-10, 10])
plt.figure("Horizontal");plt.clf()
plt.hist(horiz_change,bins=bins)
plt.xlim([0, 10])
print (depchange/3).mean(),horiz_change.mean(),calctime.mean()
