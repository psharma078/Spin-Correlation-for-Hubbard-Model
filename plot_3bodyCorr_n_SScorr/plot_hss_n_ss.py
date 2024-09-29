import matplotlib.pyplot as plt
import numpy as np
from bonds_n_coords import lat_bonds
from bonds_n_coords import coordinates
import sys
import csv

Nx=30
Ny=6
N=Nx*Ny

##hSS correlation
fname_hszsz_corr = sys.argv[1]
##Spin-Spin correlation SoSi
fname = sys.argv[2]

ij, szsz = [],[]
with open(fname_hszsz_corr, mode='r') as file:
    csv_reader = csv.reader(file)
    for row in csv_reader:
        ii = int(row[0])
        jj = int(row[1])
        ij.append([jj,ii])
        szsz.append(float(row[2]))

ij = np.array(ij)
szsz = np.array(szsz)

#for i in range(len(ij)):
#    print(ij[i], szsz[i])

fig, ax = plt.subplots(figsize=(12,4))
bond,bonds = lat_bonds(Nx,Ny)
Rxy = coordinates(Nx,Ny)

##Plot lattice geometry
def plot_lat(figo,axo):
    for k,bnd in enumerate(bonds):
        s1=bnd[0]
        s2=bnd[1]
        axo.plot([Rxy[s1][0],Rxy[s2][0]],[Rxy[s1][1],Rxy[s2][1]],lw=1,c='g',linestyle='--')
        #axo.text(Rxy[s1][0]+0.05,Rxy[s1][1],str(s1),fontsize=10)#,weight="bold",color='green')
    axo.axis("off")
    axo.axis("equal")
    #ax.text(Rxy[N][0]+0.05,Rxy[N][1],str(N),fontsize=10)
    return figo, axo

fig, ax = plot_lat(fig,ax)
##Plot hole density (for 2 hole problem)
#for i in range(1,N+1):
#    if (i!=int(N/2)+1): ax.scatter(Rxy[i][0], Rxy[i][1],s=nhole[i-1]*8000,marker='o',color='green')

##Plot hole-spin-spin correlation
bonds = np.array(bonds)
comparison = (ij[:, None, :] == bonds).all(axis=2)
indices = np.where(comparison.any(axis=1))[0]
ij = ij[indices]
szsz = szsz[indices]

scale = 30
for k,bnd in enumerate(ij):
    s1=bnd[0]
    s2=bnd[1]
    print(s1, s2, szsz[k])
    if szsz[k]<0:
        ax.plot([Rxy[s1][0],Rxy[s2][0]],[Rxy[s1][1],Rxy[s2][1]],lw=abs(szsz[k])*scale,c='red',alpha=0.7)
    else:
        ax.plot([Rxy[s1][0],Rxy[s2][0]],[Rxy[s1][1],Rxy[s2][1]],lw=abs(szsz[k])*scale,c='blue',alpha=0.7)


ij, sxy,szsz = [],[],[]
with open(fname, mode='r') as file:
    csv_reader = csv.reader(file)
    for row in csv_reader:
        ii = int(row[0])
        jj = int(row[1])
        ij.append([jj,ii])
        spm = float(row[2])
        smp = float(row[3])
        sxy.append((spm+smp)/2.)
        szsz.append(float(row[4]))

ij = np.array(ij)
sxy = np.array(sxy)
szsz = np.array(szsz)

ref_site = int(N/2)+4 #change reference site as needed
index = list((ij[:,1]==ref_site).nonzero()[0][:-1])
index.extend((ij[:,0]==ref_site).nonzero()[0][1:])
ij = np.array(ij[index])
sxy = sxy[index]
szsz = szsz[index]
ss_corr = sxy+szsz

for i in range(len(ij)):
    print(ij[i], sxy[i], szsz[i], ss_corr[i])

sites = list(np.arange(1,ref_site))
sites.extend(list(np.arange(ref_site+1,N+1)))

coord = np.array([Rxy[j] for j in sites])
x = coord[:,0]
y = coord[:,1]
u = np.zeros_like(ss_corr)
v = 7.5*np.sign(ss_corr) * np.abs(ss_corr)
#print(v)
theta = np.radians(6)

# Rotation matrix
rotation_matrix = np.array([[np.cos(theta), -np.sin(theta)],
                            [np.sin(theta),  np.cos(theta)]])

# Apply rotation to the (u, v) components
u_rotated = u * np.cos(theta) - v * np.sin(theta)
v_rotated = u * np.sin(theta) + v * np.cos(theta)
x = x - u_rotated / 2
y = y - v_rotated / 2

ax.quiver(x, y, u_rotated, v_rotated, angles='xy', scale_units='xy', scale=1, color=['b' if val > 0 else 'r' for val in ss_corr],width=0.004)

fig.savefig('hss_n_ss_corr_infU_tp_0p24.pdf',bbox_inches='tight',dpi=200)
