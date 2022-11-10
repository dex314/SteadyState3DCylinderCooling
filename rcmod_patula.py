# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 15:36:37 2016

@author: sch5926
"""

import warnings
import sys
#import numpy as np
from numpy import *
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
#import seaborn as sns

#global variables
R = 0.3
q0 = 13.7E5
h = 3410
Tc = 20
phi = 10*pi/180
psi = 90*pi/180
Q = q0*phi/h
v = 15.2
w = v/R
k = 0.03*h*R
a = v*R/(100**2)
thd = range(0,360,2)

    
ratio = [1.0, 0.999, 0.995, 0.99, 0.95, 0.9, 0.75, 0.5, 0]
rL = len(ratio)

def main():
       
    CLT = [] #define an empty set for CLT array
        
    r = np.array(ratio)*R    
    thr = np.array(thd)*pi/180

    for an in range(0,360,5):
        
        alp = an*pi/180
        gam = phi + alp
        zet = phi + alp + psi
        
        
        A0 = 1/psi
        Tsum = np.zeros((180,rL))
    
        for n in range(1,40):
            delta = sqrt(v*R*a)*(k/h/R)        
            
            x = sqrt( 1j*n*w/a )        
            Jr = exp(x*(r-R))
            JR = exp(x*(R-R))      
            
            #x = sqrt(n*w/2/a)
            #Jr = (1-1j)*exp(x*(1+1j)*(r-R))/2
            #JR = (1-1j)*exp(x*(1+1j)*(R-R))/2
            
            eth = np.exp(1j*(thr*n))
        #==============================================================    
            #3.17.16 - now this is very interesting bcs matlab reqs the -1j inorder 
            #for it to plot the function in the right direction which means there is
            #some difference that were not aware of!
        #=============================================================
                                    
            con = h/k/pi/x
            
            #print Jr.shape
            #print eth.shape
            
            eF = ( 1j/n*( np.exp(-1j*n*phi) - 1) )
            eP = ( 1j/n*( np.exp(-1j*n*zet) - np.exp(-1j*n*gam) ) )
            
            An = (eF/phi - eP/psi) * con
            A0inc = An * JR * (-1j/n*(np.exp(1j*n*zet) - np.exp(1j*n*gam)))
            A0 -= 1/psi*A0inc
            
            inc = An*Jr         
            mat = np.matmul(np.reshape(eth,(180,1)),np.reshape(inc,(1,rL)))
            
            Tsum += mat.real

        #END N LOOP-------------------------------------------------
                           
        CLT.append( A0 )
        TF = Tsum + A0
        if an == 0:
            Tf0 = np.reshape(TF,(180,rL))
        elif an == 15:
            Tf15 = np.reshape(TF,(180,rL))
        elif an == 30:
            Tf30 = np.reshape(TF,(180,rL))
        elif an == 45:
            Tf45 = np.reshape(TF,(180,rL))
            
        #if an == 0 or an == 15 or an == 30 or an == 45:
            #fign = an/15 + 1
            #plt.figure(fign)
            #plt.plot(np.reshape(thd,(180,1)),np.reshape(TF,(180,rL)))

    #END AN LOOP----------------------------------------------------
        
    #Plotting Section-----------------------------------------------    
    f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex = True, sharey = True)
    ax1.plot(np.reshape(thd,(180,1)),Tf0), ax1.set_title('Moving Alpha Angle')
    ax2.plot(np.reshape(thd,(180,1)),Tf15)
    ax3.plot(np.reshape(thd,(180,1)),Tf30)
    ax4.plot(np.reshape(thd,(180,1)),Tf45)
    f.subplots_adjust(hspace=0) #decrease plot spacing
    plt.setp([ax.get_xticklabels() for ax in f.axes[:-1]], visible=False) #minimizes axis ticks

    
    plt.figure(2)
    plt.plot(CLT), plt.title('Center Line Temperature Over Alpha Change')    
    ylim([0, 1])
    
    rm, tm = meshgrid(r,thd)
    #plt.figure(3)
    #hcont = plt.contourf(rm,tm,Tf45)   

    sfig = plt.figure(4)    
    ax2 = sfig.gca(projection='3d')
    surf = ax2.plot_surface(rm,tm,Tf45, rstride=1, cstride=1, cmap=cm.jet,
                       linewidth=0, antialiased=False)
    
if __name__ == '__main__':
    main()              
        
        
