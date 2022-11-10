# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 09:44:16 2016

@author: sch5926
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 15:36:37 2016

@author: sch5926
"""

import warnings
import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

def pol2cart(rho, phi, Nr, Np):
    x = np.reshape(rho,(Nr,1)) * np.reshape(np.cos(phi),(1,Np))
    y = np.reshape(rho,(Nr,1)) * np.reshape(np.sin(phi),(1,Np))
    return(x, y)

#global variables
R = 0.3
q0 = 13.7E5
h = 3410
Tc = 20
phi = 10*pi/180
#psi = 90*pi/180
psi = np.array([90, 0])*pi/180
Q = q0*phi/h
v = 1.0
w = v/R
k = 0.03*h*R
a = 4.6E-4 #v*R/(100**2)
thd = range(0,362,2)
thL = len(thd)

    
ratio = [1.0, 0.999, 0.998, 0.997, 0.996, 0.995, 0.99, 0.95, 0.9, 
         0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0]
rL = len(ratio)

rlen = 2.032
L = rlen/2
sw = L*0.75
z0 = np.linspace(0,L,29)
z = np.array(z0)
zL = len(z)
hL = np.array([1*L, 1*L])

def main():
       
    CLT = [] #define an empty set for CLT array
        
    r = np.array(ratio)*R    
    thr = np.array(thd)*pi/180

    for an in range(0,360,15):
         
        alp1 = an*pi/180
        alp2 = 0*pi/180
        G1 = phi + alp1
        Z1 = phi + alp1 + psi[0]
        G2 = phi + alp1 + psi[0] + alp2
        Z2 = phi + alp1 + psi[0] + alp2 + psi[1]        
            
        A0 = 0
        A0m = np.zeros((rL,zL))
        Tsum = np.zeros((rL*zL,thL))
    
        con = 2*h/k/pi/L
        sumPL = sum(psi * hL)
        for m in range(0,2*(zL-1)):
            
            b = (2*m+1)*pi/2/L 
            ZS = sp.sin(b*sw)/b
            Zm = sp.cos(b*z)
            Jz = sp.i0(b*r)
            JZ = sp.i0(b*R)
            
            for n in range(1,40):
                delta = sqrt(v*R*a)*(k/h/R)        
                
                lam = 1j*n*w/a
                x = sqrt(lam + b**2)
                
                Jr = np.exp(x*r)/np.sqrt(2*pi*x)
                JR = np.exp(x*R)/np.sqrt(2*pi*x)      
                            
                eth = np.exp(1j*(thr*n))
                #==============================================================    
                #3.17.16 - now this is very interesting bcs matlab reqs the -1j inorder 
                #for it to plot the function in the right direction which means there is
                #some difference that were not aware of!
                #=============================================================
                
                eF = ( 1j/n*( np.exp(-1j*n*phi) - 1) )
                
                eP1 = ( 1j/n*( np.exp(-1j*n*Z1) - np.exp(-1j*n*G1) ) )                
                eP2 = ( 1j/n*( np.exp(-1j*n*Z2) - np.exp(-1j*n*G2) ) )  
                ePsum = (hL[0]*eP1 + hL[1]*eP2)
                
                ePn1 = (-1j/n*( np.exp(1j*n*Z1) - np.exp(1j*n*G1) )) 
                ePn2 = (-1j/n*( np.exp(1j*n*Z2) - np.exp(1j*n*G2) )) 
                ePnsum = (hL[0]*ePn1 + hL[1]*ePn2)
                
                
                An = (eF/phi - ePsum/sumPL) * ZS * con/(x*JR)
                A0inc = (An * JR * ePnsum) / sumPL
                A0 -= A0inc.real
                
                incn = An*Jr         
                mat_rz = np.dot(np.reshape(incn,(rL,1)),
                                   np.reshape(Zm,(1,zL)))
                                
                mat = np.dot( np.reshape(mat_rz,(rL*zL,1)), 
                             np.reshape(eth,(1,thL)) )
                
                Tsum += mat.real
    
            #END N LOOP-------------------------------------------------
            incm = Jz/JZ*ZS
            A0minc = 2 * (1/sumPL+A0) * np.dot(np.reshape(incm,(rL,1)),
                          np.reshape(Zm,(1,zL)))
            A0m += A0minc.real
            
        #END M LOOP-----------------------------------------------------
        CLT.append( A0m[0,0] )
        
        A0mf = np.dot(np.reshape(A0m,(rL*zL,1)), np.ones((1,thL)))
        
        TF = Tsum + A0mf
    
        t_f = np.reshape(TF,(rL,zL,thL))                   
        #print t_f[0,0,4]
        
        t_rt = t_f[:,0,:]
        t_rz = np.array(t_f[:,:,4])
        
        if an == 0:
            Tfrt0 = t_rt
            Tfrz0 = t_rz
        elif an == 15:
            Tfrt15 = t_rt
            Tfrz15 = t_rz
        elif an == 30:
            Tfrt30 = t_rt
            Tfrz30 = t_rz
        elif an == 45:
            Tfrt45 = t_rt
            Tfrz45 = t_rz
           
#        plt.figure(1)
#        plt.plot(np.transpose(thd),np.transpose(t_rt))
#        plt.title('Thermal Distribution in R and Theta')
#        
#        plt.figure(2)
#        plt.plot(np.transpose(z),np.transpose(t_rz))
#        plt.title('Thermal Distribution in R and Z')
        
    #END an ALPHA LOOP--------------------------------------------------
 
        
    #%%Plotting Section-----------------------------------------------    
    frt, (ar1, ar2, ar3, ar4) = plt.subplots(4, sharex = True, sharey = True)
    ar1.plot(np.transpose(thd),np.transpose(Tfrt0)), ar1.set_title('Moving Alpha Angle - R/Theta')
    ar2.plot(np.transpose(thd),np.transpose(Tfrt15))
    ar3.plot(np.transpose(thd),np.transpose(Tfrt30))
    ar4.plot(np.transpose(thd),np.transpose(Tfrt45))
    frt.subplots_adjust(hspace=0) #decrease plot spacing
    plt.setp([ar.get_xticklabels() for ar in frt.axes[:-1]], visible=False) #minimizes axis ticks
    
    frz, ((az1, az2), (az3, az4)) = plt.subplots(2, 2, sharex = 'col', sharey = 'row')
    az1.plot((z),np.transpose(Tfrz0)), az1.set_title('Moving Alpha Angle - R/Z')
    az2.plot((z),np.transpose(Tfrz15))
    az3.plot((z),np.transpose(Tfrz30))
    az4.plot((z),np.transpose(Tfrz45))
    
    plt.figure()
    plt.plot(CLT), plt.title('Center Line Temperature Over Alpha Change')    
    #ylim([0, 1])
    
    xc,yc = pol2cart(r,thr,rL,thL)
    
    rm, tm = np.meshgrid(r,thd)
    print rm.shape
    sfig = plt.figure()    
    ax2 = sfig.gca(projection='3d')
    
#    surf = ax2.plot_surface(rm,tm,np.transpose(Tfrt45), rstride=1, cstride=1, cmap=cm.jet,
#                       linewidth=0, antialiased=False)
    surf = ax2.plot_surface(xc,yc,Tfrt45, rstride=1, cstride=1, cmap=cm.jet,
                       linewidth=0, antialiased=False)
    
    #Xc,Yc = np.meshgrid(xc,yc)
    #print xc.shape
    #print yc.shape    
    plt.figure()
    plt.subplot(2,2,1)
    plt.pcolor(xc,yc,Tfrt0, cmap=cm.jet)
    plt.colorbar()
    plt.subplot(2,2,2)
    plt.pcolor(xc,yc,Tfrt15, cmap=cm.jet)
    plt.colorbar()
    plt.subplot(2,2,3)
    plt.pcolor(xc,yc,Tfrt30, cmap=cm.jet)
    plt.colorbar()
    plt.subplot(2,2,4)
    plt.pcolor(xc,yc,Tfrt45, cmap=cm.jet)
    plt.colorbar()
    #surf = ap1.plot_surface(Xc, Yc, Tfrt0, rstride=1, cstride=1, cmap=cm.jet,
    #                   linewidth=0, antialiased=False)
   

#    
#    rm, tm = meshgrid(r,thd)
#    #plt.figure(3)
#    #hcont = plt.contourf(rm,tm,Tf45)   
#
#    sfig = plt.figure(4)    
#    ax2 = sfig.gca(projection='3d')
#    surf = ax2.plot_surface(rm,tm,Tf45, rstride=1, cstride=1, cmap=cm.jet,
#                       linewidth=0, antialiased=False)
    
#%%END
if __name__ == '__main__':
    main()              
    

        
        
