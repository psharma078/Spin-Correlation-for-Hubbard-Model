import matplotlib.pyplot as plt
import numpy as np

def lat_bonds(Nx,Ny):
    N = Nx*Ny
    bond=[]
    bonds = []
    for i in range(1,N):
        if i%Ny==1 and i<N-Ny:
            b=[i,i+1]
            bond.append(b)
            bonds.append(b)
            b=[i,i+Ny]
            bonds.append(b)
            bond.append(b)
            b=[i,i+Ny+1]
            bonds.append(b)
            bond.append(b)
            b=[i,i+Ny-1]
            #bonds.append(b)
            bond.append(b)
            b=[i,i+2*Ny-1]
            #bonds.append(b)
            bond.append(b)
   
        elif (i%Ny==0 and i<N):
            b=[i,i+Ny]
            bond.append(b)
            bonds.append(b)

        elif i>N-Ny:
            b=[i,i+1]
            bond.append(b)
            bonds.append(b)
            if i==N-Ny+1:
                b=[i,i+Ny-1]
                bond.append(b)

        elif i%2==1:
            b=[i,i+1]
            bond.append(b)
            bonds.append(b)
            b=[i,i+Ny-1]
            bond.append(b)
            bonds.append(b)
            b=[i,i+Ny]
            bond.append(b)
            bonds.append(b)
            b=[i,i+Ny+1]
            bond.append(b)
            bonds.append(b)

        elif i%2==0:
            b=[i,i+1]
            bond.append(b)
            bonds.append(b)
            b=[i,i+Ny]
            bond.append(b)
            bonds.append(b)
    return bond, bonds     #bond for cylindrical bc and bonds for open

#for i  in range(len(bond)):
#    print(bond[i])

#finding coordinates
def coordinates(Nx,Ny):
    N=Nx*Ny
    a1 = np.array([1,0])
    a2 = np.array([0.5,np.sqrt(3)/2.0])

    Rxy = {}

    for j in range(Ny,0,-1):
        for i in range(1,Nx+1):
            label = j+(i-1)*Ny
            if j == 6:
                Rxy[label] = (i-1)*a1
            elif j==5:
                Rxy[label] = (i-1)*a1 + a2
            elif j==4:
                Rxy[label] = (i-2)*a1 + 2*a2
            elif j==3:
                Rxy[label] = (i-2)*a1 + 3*a2
            elif j==2:
                Rxy[label] = (i-3)*a1 + 4*a2
            else:
                Rxy[label] = (i-3)*a1 + 5*a2
    return Rxy

'''
bond,bonds = lat_bonds()
Rxy = coordinates()
for k,bnd in enumerate(bonds):
    s1=bnd[0]
    s2=bnd[1]
    
    plt.plot([Rxy[s1][0],Rxy[s2][0]],[Rxy[s1][1],Rxy[s2][1]],lw=1,c='k')
    plt.text(Rxy[s1][0]+0.05,Rxy[s1][1],str(s1),fontsize=10)#,weight="bold",color='green')
plt.text(Rxy[36][0]+0.05,Rxy[36][1],str(36),fontsize=10)
plt.axis("off")
#plt.savefig('lattice_XC6.pdf',dpi=600,bbox_inches='tight')
#print(bonds)


plt.show()
'''
