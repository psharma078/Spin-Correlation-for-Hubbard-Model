from bonds_n_coords import lat_bonds
from bonds_n_coords import coordinates
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import matplotlib.pyplot as plt
import sys
import os
import numpy as np


fname = sys.argv[1]
Nx = 6
Ny = 6
N = Nx*Ny
lineNum = int(N*(N+1)/2)
command = "grep -A"+str(lineNum)+" \"i j <Si^z Sj^z>\" "+fname+" > tmp"
os.system(command)

SiSj = [[],[],[],[],[]]
f=open("tmp",'r')
for l,line in enumerate(f):
    if (l>=1):
        line=line.strip()
        line=line.split(" ")
        #print(line[0], line[1], line[2], line[3])
        SiSj[0].append(int(line[0]))
        SiSj[1].append(int(line[1]))
        SiSj[2].append(float(line[2]))
        SiSj[3].append(float(line[3]))
        SiSj[4].append(float(line[4]))
f.close()

S2 = 0
for i in range(len(SiSj[0])):
    m = SiSj[0][i]
    n = SiSj[1][i]
    szz = SiSj[2][i]
    spm = SiSj[3][i]
    smp = SiSj[4][i]
    if m==n:
        S2 +=  szz + 0.5*(spm+smp)
    else: S2 += 2*(szz + 0.5*(spm+smp))

print("S^2= ",S2)
command = "grep \"Energy using overlap\" "+fname
energy = os.system(command)
print(energy)

#------------------------------------
bonds = lat_bonds(Nx,Ny)
Rxy = coordinates(Nx,Ny)

KX = []
KY = []
for n in range(-Nx,Nx+1):
    for m in range(-Ny,Ny+1):
        Kx = 1.0/Nx*n       #in the unit of 2pi
        Ky = (2.0/np.sqrt(3))*m/Ny
        KX.append(Kx)
        KY.append(Ky)

def str_fac(kx,ky,r,corr):
    SFxx = 0
    SFzz = 0
    for i in range(len(corr[0])):
        index_i = corr[0][i]
        index_j = corr[1][i]
        szsz = corr[2][i]
        sxsx = 1.0/4.0*(corr[3][i]+corr[4][i])
        
        r_ij = r[index_j]-r[index_i]
        #kr_ij = np.dot(k_vec,r_ij)
        kr_ij = kx*r_ij[0]+ky*r_ij[1]
        if (index_i==index_j):
            SFzz += szsz
            SFxx += sxsx
        else:
            SFzz += 2*np.cos(kr_ij*2*np.pi)*szsz
            SFxx += 2*np.cos(kr_ij*2*np.pi)*sxsx

    return SFzz/N, SFxx/N

SFz, SFx = str_fac(0,0,Rxy,SiSj)
print("S(q=0)=", SFz, SFx, (SFz+2*SFx))

#kx,ky = np.meshgrid(KX,KY)

kx = np.reshape(KX,(2*Nx+1,2*Ny+1))
ky = np.reshape(KY,(2*Nx+1,2*Ny+1))

Sfzz,Sfxx = str_fac(kx,ky,Rxy,SiSj)
print(Sfzz.min(), Sfzz.max())
print(Sfxx.min(), Sfxx.max())
#stot = Sfzz+Sfxx
#print(stot.min(), stot.max())
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#TEXT SETTING
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from matplotlib import rc
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 8}
rc('font', **font)
rc('text', usetex=True)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#PLOT STATIC STRUCTURE FACTOR OVER THE K SPACE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fig,(ax1,ax2) = plt.subplots(nrows=1, ncols=2,sharey=True,subplot_kw={'xticks': [-1,-0.5,0,0.5,1.0], 'yticks': [-1,-0.5,0,0.5,1.0]},figsize=(4.2,2))
fig.subplots_adjust(wspace=0.17)
im1 = ax1.contourf(kx,ky,Sfzz,300,extent=[-1, 1, -1, 1],cmap=plt.cm.gist_heat)
ax1_divider = make_axes_locatable(ax1)
for im in im1.collections:
    im.set_edgecolor("face")
# add an axes to the right of the main axes.
cax1 = ax1_divider.append_axes("top", size="5%", pad="2%")

length=Sfzz.max()-Sfzz.min()
#tick1=[round(Sfzz.min()+x*(length)*1.0/4,2) for x in range(5)]
tick1=[round(Sfzz.min(),2),round((Sfzz.min()+Sfzz.max())/2,2),round(Sfzz.max(),2)-0.002]

cb1 = fig.colorbar(im1, cax=cax1,orientation="horizontal",ticks=tick1)
cax1.xaxis.set_ticks_position("top")
ax1.set_title(r'$S^{zz}(\mathbf q)$',y=1.19)
ax1.set_ylim([-1,1])
ax1.set_xlabel(r'$q_x/2\pi$',fontsize = 12)
ax1.set_ylabel(r'$q_y/2\pi$',rotation=90, position=(-1,0.5),fontsize = 12)

#-----------------------------------------------------------
im2 = ax2.contourf(kx,ky,Sfxx,300,extent=[-1, 1, -1, 1],cmap=plt.cm.gist_heat)
ax2_divider = make_axes_locatable(ax2)

for im in im2.collections:
    im.set_edgecolor("face")
# add an axes above the main axes.
cax2 = ax2_divider.append_axes("top", size="5%", pad="2%")
#tick2=[round(x*(Sfxx.max()-0.005)*1.0/5,2) for x in range(6)]
tick2=[round(Sfxx.min()+0.01,2),round((Sfxx.min()+Sfxx.max())/2,2),round(Sfxx.max(),2)-0.01]
cb2 = fig.colorbar(im2, cax=cax2, orientation="horizontal",ticks=tick2)
cax2.xaxis.set_ticks_position("top")
ax2.set_title(r'$S^{xx}(\mathbf q)$',y=1.19)
ax2.set_xlabel(r'$q_x/2\pi$',fontsize = 12)

#corners of hexagonal BZ
coord=[[2.0/3,0],[1.0/3,-1.0/np.sqrt(3)],[-1.0/3,-1.0/np.sqrt(3)],[-2.0/3,0],[-1.0/3,1.0/np.sqrt(3)],[1.0/3,1.0/np.sqrt(3)]]
bonds=[[1,2],[2,3],[3,4],[4,5],[5,6],[6,1]]

for j,bond in enumerate(bonds):
    site1=bond[0]-1
    site2=bond[1]-1
    ax1.plot([coord[site1][0],coord[site2][0]],[coord[site1][1],coord[site2][1]],lw=0.9,c='y',linestyle='--')
    ax2.plot([coord[site1][0],coord[site2][0]],[coord[site1][1],coord[site2][1]],lw=0.9,c='y',linestyle='--')

#plt.savefig("XC6_6x6_Bragg_peak_nup=29_ndn=29.pdf",dpi=600,bbox_inches='tight')
plt.show()

