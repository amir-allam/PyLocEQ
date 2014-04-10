from matplotlib import pyplot as plt
from numpy import arange

def plot_z_slice(data,qx=None,qy=None,qz=None,iz=0,fignum=None,if_faults=False,cscale=None):
    #Plot a map-view slice through the 3D matrix data at index iz
    # qx,qy,qz are the coordinates along the dimensions of the matrix data
    if qx is None:
        qx=arange(0,data.shape[1])
    if qy is None:
        qy=arange(0,data.shape[0])
    if qz is None:
        qz=arange(0,data.shape[2])
    figid=plt.figure(fignum)
    figid.clf();
    if if_faults:
        from misc_tools import load_faults
        faultlon,faultlat=load_faults()
        plt.plot(faultlon,faultlat,'k')
    plt.axis([min(qx),max(qx),min(qy),max(qy)])
    pc_id=plt.pcolor(qx,qy,data[:,:,iz] )
    if cscale is not None:
        plt.clim(cscale)
    figid.colorbar(pc_id)
    figid.show()
