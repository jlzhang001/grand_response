import scipy.io as scio
import numpy as np 
import cmath
import os

path = './EFL_all.mat'
data = scio.loadmat(path)
print(type(data))
print(data.keys())
Theta=data['Theta']
Phi=data['Phi']
freq_All=data['freq_ALL']
EffL_all=data['EffL_all']
print(np.shape(EffL_all))

#nangles=int(91*(355/5+1))
#nangles=int(181*(355/5+1))
nangles=int(181*361)
realimp=np.zeros((freq_All.shape[1]))
reactance=np.zeros((freq_All.shape[1]))
theta=np.zeros((len(realimp),nangles))
phi=np.zeros((len(realimp),nangles))

realimpbis=np.zeros((len(realimp),nangles))
freqbis=np.zeros((len(realimp),nangles))
reactancebis=np.zeros((len(realimp),nangles))
EffLx=np.zeros((3,len(realimp),nangles),dtype=complex)
EffLy=np.zeros((3,len(realimp),nangles),dtype=complex)
EffLz=np.zeros((3,len(realimp),nangles),dtype=complex)
EffL_r=np.zeros((3,len(realimp),nangles),dtype=complex)
EffL_theta=np.zeros((3,len(realimp),nangles),dtype=complex)
EffL_phi=np.zeros((3,len(realimp),nangles),dtype=complex)
phase_theta=np.zeros((3,len(realimp),nangles))
phase_phi=np.zeros((3,len(realimp),nangles))

print("for loop: ")
for i in range(0,freq_All.shape[1]):
    for j in range(0,Phi.shape[1]):
        #for k in range(0,Theta.shape[1]):
            for nant in range(0,3):
                #http://dynref.engr.illinois.edu/rvs.html note, theta is azimuth, phi is zenith!!!
                #formula: cartesian to spherical
                #Er = Ex*sin(theta)*cos(phi)+ Ey*sin(theta)*sin(phi) + Ez*cos(Theta)
                #Etheta=Ex*cos(theta)*cos(phi)+ Ey*cos(theta)*sin(phi) - Ez*sin(Theta)
                #Ephi= -Ex*sin(phi)+Ey*cos(phi))
                EffL_r[nant,i,j*Theta.shape[1]:j*Theta.shape[1]+Theta.shape[1]]=EffL_all[i,nant,0,:,j]*np.sin(Theta[0,:]*np.pi/180.)*np.cos(Phi[0,j]*np.pi/180.) + EffL_all[i,nant,1,:,j]*np.sin(Theta[0,:]*np.pi/180.)*np.sin(Phi[0,j]*np.pi/180.)+EffL_all[i,nant,2,:,j]*np.cos(Theta[0,:]*np.pi/180.)
                EffL_theta[nant,i,j*Theta.shape[1]:j*Theta.shape[1]+Theta.shape[1]]=EffL_all[i,nant,0,:,j]*np.cos(Theta[0,:]*np.pi/180.)*np.cos(Phi[0,j]*np.pi/180.) + EffL_all[i,nant,1,:,j]*np.cos(Theta[0,:]*np.pi/180.)*np.sin(Phi[0,j]*np.pi/180.)-EffL_all[i,nant,2,:,j]*np.sin(Theta[0,:]*np.pi/180.)
                phase_theta[nant,i,j*Theta.shape[1]:j*Theta.shape[1]+Theta.shape[1]]=np.degrees(np.angle(EffL_theta[nant,i,j*Theta.shape[1]:j*Theta.shape[1]+Theta.shape[1]]))

                EffL_phi[nant,i,j*Theta.shape[1]:j*Theta.shape[1]+Theta.shape[1]]=EffL_all[i,nant,0,:,j]*(-np.sin(Phi[0,j]*np.pi/180.))+EffL_all[i,nant,1,:,j]*np.cos(Phi[0,j]*np.pi/180.)
                phase_phi[nant,i,j*Theta.shape[1]:j*Theta.shape[1]+Theta.shape[1]]=np.degrees(np.angle(EffL_phi[nant,i,j*Theta.shape[1]:j*Theta.shape[1]+Theta.shape[1]]))

for j in range(0,nangles):
    freqbis[:,j]=freq_All*1000.
    realimpbis[:,j]=realimp
    reactancebis[:,j]=reactance

for i in range(0,freq_All.shape[1]):
    for j in range(0,Phi.shape[1]):
        theta[i,j*Theta.shape[1]:(j+1)*Theta.shape[1]]=Theta
        phi[i,j*Theta.shape[1]:(j+1)*Theta.shape[1]]=Phi[0,j]


#np.save('./antenna_leff_x', (freqbis,realimpbis,reactancebis,theta,phi,lefftheta,leffphi,phasetheta,phasephi))
np.save('./antenna_leff_x2', (freqbis,realimpbis,reactancebis,theta,phi,abs(EffL_theta[0,:,:]),abs(EffL_phi[0,:,:]),phase_theta[0,:,:],phase_phi[0,:,:]))

np.save('./antenna_leff_y2', (freqbis,realimpbis,reactancebis,theta,phi,abs(EffL_theta[1,:,:]),abs(EffL_phi[1,:,:]),phase_theta[1,:,:],phase_phi[1,:,:]))

np.save('./antenna_leff_z2', (freqbis,realimpbis,reactancebis,theta,phi,abs(EffL_theta[2,:,:]),abs(EffL_phi[2,:,:]),phase_theta[2,:,:],phase_phi[2,:,:]))

print("saved!")

#plot=0 #no plot
plot=1 #plot

if plot<=0:
   sys.exit()

import matplotlib.pyplot as plt

myEffL_all=np.zeros((221,3,3,181,361),dtype=complex)

for i in range(0,freq_All.shape[1]):
    for j in range(0,Phi.shape[1]):
        #for k in range(0,Theta.shape[1]):
            for nant in range(0,3):
                #http://dynref.engr.illinois.edu/rvs.html#rvs-et-d
                #formula: spherical to cartesian
                #Ex = Er*sin(theta)*cos(phi)- Etheta*cos(theta)*cos(phi) - Ephi*sin(phi)
                #Ey = Er*sin(theta)*sin(phi)- Etheta*cos(theta)*sin(phi) + Ephi*cos(phi) 
                #Ez = Er*cos(theta)-Etheta*sin(theta))
                 
                myEffL_all[i,nant,0,:,j]=EffL_r[nant,i,j*Theta.shape[1]:j*Theta.shape[1]+Theta.shape[1]]*np.cos(Phi[0,j]*np.pi/180.)*np.sin(Theta[0,:]*np.pi/180.) - EffL_phi[nant,i,j*Theta.shape[1]:j*Theta.shape[1]+Theta.shape[1]]*np.sin(Phi[0,j]*np.pi/180.) + np.cos(Theta[0,:]*np.pi/180.)*np.cos(Phi[0,j]*np.pi/180.)*EffL_theta[nant,i,j*Theta.shape[1]:j*Theta.shape[1]+Theta.shape[1]]
                myEffL_all[i,nant,1,:,j]=EffL_r[nant,i,j*Theta.shape[1]:j*Theta.shape[1]+Theta.shape[1]]*np.sin(Phi[0,j]*np.pi/180.)*np.sin(Theta[0,:]*np.pi/180.) + EffL_phi[nant,i,j*Theta.shape[1]:j*Theta.shape[1]+Theta.shape[1]]*np.cos(Phi[0,j]*np.pi/180.) + np.cos(Theta[0,:]*np.pi/180.)*np.sin(Phi[0,j]*np.pi/180.)*EffL_theta[nant,i,j*Theta.shape[1]:j*Theta.shape[1]+Theta.shape[1]]
                myEffL_all[i,nant,2,:,j]=EffL_r[nant,i,j*Theta.shape[1]:j*Theta.shape[1]+Theta.shape[1]]*np.cos(Theta[0,:]*np.pi/180.)-np.sin(Theta[0,:]*np.pi/180.)*EffL_theta[nant,i,j*Theta.shape[1]:j*Theta.shape[1]+Theta.shape[1]]


thetabin=45
phibin=45
freqbin=40

thetaValue=Theta[0][thetabin]
phiValue=Phi[0][phibin]
freqValue=freq_All[0][freqbin]

print("Plotting for Freq",freqValue*1000,"Theta",thetaValue,"Phi",phiValue)

def PhaseSpacePlot(fig,ax,PlotTitle,X,Y,Z):
    tmp=ax.pcolormesh(X, Y, Z.T,shading="auto",cmap="bwr",vmax=1.5,vmin=-1.5)
    tmp=fig.colorbar(tmp, ax=ax)
    #tmp=ax.set_aspect('equal')
    tmp=ax.set_title(PlotTitle)
    tmp=ax.set_ylabel('$\\Phi [deg]$')
    tmp=ax.set_xlabel('$\\Theta [deg]$')
    #for secondary axis
    #lambdaZenith = lambda LSZ: np.rad2deg(np.arccos(np.power(10,-LSZ)))
    #lambdaLSZ = lambda Zenith: np.log10(1.0/np.cos(np.deg2rad(Zenith)))
    #ax2 = ax.secondary_yaxis('right',functions=(lambdaZenith,lambdaLSZ))
    #ax2.set_yticks(Zeniths[0:])
    #ax2.set_yticklabels([str(round(float(label), 1)) for label in Zeniths[0:]]) 
    #tmp=ax2.set_ylabel("$Zenith [deg]$")
    #
fig7 = plt.figure(72,figsize=(7,5), facecolor='w', edgecolor='k')

magicbin=0

ax1=fig7.add_subplot(331)
yvalues=Phi[0]
xvalues=Theta[0][0:90]
zvalues=myEffL_all[freqbin,magicbin,0,0:90,:].real
name = 'Response at '+str(freqValue*1000)+" Mhz, Magic "+str(magicbin)+" X"
PhaseSpacePlot(fig7,ax1,name,xvalues,yvalues,zvalues)

ax1=fig7.add_subplot(332)
yvalues=Phi[0]
xvalues=Theta[0][0:90]
zvalues=myEffL_all[freqbin,magicbin,1,0:90,:].real
name = 'Response at '+str(freqValue*1000)+" Mhz, Magic "+str(magicbin)+" Y"
PhaseSpacePlot(fig7,ax1,name,xvalues,yvalues,zvalues)

ax1=fig7.add_subplot(333)
yvalues=Phi[0]
xvalues=Theta[0][0:90]
zvalues=myEffL_all[freqbin,magicbin,2,0:90,:].real
name = 'Response at '+str(freqValue*1000)+" Mhz, Magic "+str(magicbin)+" Z"
PhaseSpacePlot(fig7,ax1,name,xvalues,yvalues,zvalues)

magicbin=1
ax1=fig7.add_subplot(334)
yvalues=Phi[0]
xvalues=Theta[0][0:90]
zvalues=myEffL_all[freqbin,magicbin,0,0:90,:].real
name = 'Response at '+str(freqValue*1000)+" Mhz, Magic "+str(magicbin)+" X"
PhaseSpacePlot(fig7,ax1,name,xvalues,yvalues,zvalues)

magicbin=1
ax1=fig7.add_subplot(334)
yvalues=Phi[0]
xvalues=Theta[0][0:90]
zvalues=myEffL_all[freqbin,magicbin,0,0:90,:].real
name = 'Response at '+str(freqValue*1000)+" Mhz, Magic "+str(magicbin)+" X"
PhaseSpacePlot(fig7,ax1,name,xvalues,yvalues,zvalues)

ax1=fig7.add_subplot(335)
yvalues=Phi[0]
xvalues=Theta[0][0:90]
zvalues=myEffL_all[freqbin,magicbin,1,0:90,:].real
name = 'Response at '+str(freqValue*1000)+" Mhz, Magic "+str(magicbin)+" Y"
PhaseSpacePlot(fig7,ax1,name,xvalues,yvalues,zvalues)

ax1=fig7.add_subplot(336)
yvalues=Phi[0]
xvalues=Theta[0][0:90]
zvalues=myEffL_all[freqbin,magicbin,2,0:90,:].real
name = 'Response at '+str(freqValue*1000)+" Mhz, Magic "+str(magicbin)+" Z"
PhaseSpacePlot(fig7,ax1,name,xvalues,yvalues,zvalues)

magicbin=2
ax1=fig7.add_subplot(337)
yvalues=Phi[0]
xvalues=Theta[0][0:90]
zvalues=myEffL_all[freqbin,magicbin,0,0:90,:].real
name = 'Response at '+str(freqValue*1000)+" Mhz, Magic "+str(magicbin)+" X"
PhaseSpacePlot(fig7,ax1,name,xvalues,yvalues,zvalues)

ax1=fig7.add_subplot(338)
yvalues=Phi[0]
xvalues=Theta[0][0:90]
zvalues=myEffL_all[freqbin,magicbin,1,0:90,:].real
name = 'Response at '+str(freqValue*1000)+" Mhz, Magic "+str(magicbin)+" Y"
PhaseSpacePlot(fig7,ax1,name,xvalues,yvalues,zvalues)

ax1=fig7.add_subplot(339)
yvalues=Phi[0]
xvalues=Theta[0][0:90]
zvalues=myEffL_all[freqbin,magicbin,2,0:90,:].real
name = 'Response at '+str(freqValue*1000)+" Mhz, Magic "+str(magicbin)+" Z"
PhaseSpacePlot(fig7,ax1,name,xvalues,yvalues,zvalues)

plt.show()
 
