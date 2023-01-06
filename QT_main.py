#!/usr/bin/env python
# coding: utf-8

# In[386]:


############ Quantum tunneling, a Numerical Study
########### Afsar Reja SSCU PhD
###########  Reproducing the paper: Simulation in Quantum tunneling by Kevin Smith and Guy Blaylock
##########for the course project  PH 354


import numpy as np
from math import exp, pi, pow,sqrt
from matplotlib import pyplot as plt

##self decleared values (IN SI UNITS)

x0=-14*1e-6 #meter
k0=2*1e6   #meter inverse
sigmak = 0.7*1e6;
#sigmak=0.5 #micrometer^-1
sigmax=1/(2*sigmak)
hbar = 1.054*1e-34
m = 9.1*1e-31
E0 = hbar**2*k0**2/(2*m)
dt = 0.05*1e-9 




## conversion factors from si to natural

lencon = 1.0/1.97e-7
kcon = 1.97e-7
econ = 1.0/1.60217657e-19
tcon = 1.0/6.58e-16
mcon = 1.0/1.7826e-36

### these values are in antural units
x0 = x0*lencon
k0 = k0*kcon
sigmax = sigmax * lencon
hbar = 1.0
m = m*mcon
E0 = E0*econ
dt = dt*tcon
print('x0=',x0)
print('E0=',E0)
print("sigmax=",sigmax)
##### we will do everything in natural units

#Wave function initialization
PSI=[]
PSI_2=[]

def psi(x):
    y = np.exp(complex(0,k0*x))*np.exp(-pow((x-x0),2)/(4*sigmax**2))
    y=y/pow(2*pi*sigmax**2,0.25)
    return y
  
xwidth = 80*1e-6 ##micrometer
xwidthR = 100*1e-6 ##micrometer
xwidth = xwidth*lencon
xwidthR = xwidthR*lencon
xgrid = 0.1*1e-6
dx=xgrid
xgrid = xgrid*lencon
X=np.arange(-xwidth,xwidthR,xgrid) #micro meter
print('xwidth=',xwidth)
print('xgrid',xgrid)


for x1 in X:
    y=psi(x1)
    PSI.append(y)
    PSI_2.append(abs(y)**2)
   # print(y,"\t",abs(y)**2)

#print(PSI_2)
plt.plot(X,PSI_2,label="PSI_2 in natural space")
plt.legend()
plt.show()


PSI_c2=np.array(PSI_2)*lencon*1e-6 #micrometer inverse
Xn=np.array(X)*1e6/lencon    #micrometer
plt.plot(Xn,PSI_c2,label="PSI_2 in x space")
plt.ylabel("$||\psi(x)||^2$")
plt.xlabel("x($\mu$m)")
plt.legend()
plt.show()



###########k space  for FFT
#### k space grid points
# N=len(X)
# dp = 2*np.pi*hbar/(N*xgrid);   # k-space width of grid point
# P=np.arange(0,N)*dp-np.average(np.arange(0,N)*dp)
# P=P*lencon*1e-6


# PHI=[]
# PHI_2=[]

# #PSI_c=np.array(PSI)*(lencon*1e-6)**0.5 #for FFT

# PHI=np.fft.fft(np.array(PSI))
# PHI_2=(np.absolute(PHI))**2
# PHI_2=PHI_2/(lencon*1e-6)

# plt.xlim(left=0,right=400)
# plt.plot(PHI_2,label="PHI_2 in k space")
# plt.legend()


# In[370]:


#potential initialization
a = 1e-6  #micrometer
a = a*lencon

V=[]
for x1 in X:
    if abs(x1)>a:
        v=0*econ
    else:
        #v=0
        v=2*E0
    V.append(v)
        
plt.plot(X,V)
#plt.plot(X,PSI_2)
#plt.axhline(y=.04,color='r')


# In[371]:


### plotting potential and PSI_2 on same plot


plt.plot(Xn,PSI_c2,label="PSI_2 in x space")
plt.ylabel("$||\psi(x)||^2$")
plt.xlabel("x ( $\mu$m )")
plt.legend()
plt.xlim(-25,25)
plt.axvline(x=-a/lencon,color="r")
plt.axvline(x=a/lencon,color="r")
plt.show()


# In[372]:


################################### A and B and C matrix formulation
N=len(X)
#N=5
dx = xgrid
A=np.zeros([N,N],dtype=complex)
B=np.zeros([N,N],dtype=complex)
C=np.zeros([N,N],dtype=complex)

alpha=complex(0,(hbar*dt/(2*m*dx*dx)))
            
for i in range(N):
    if i==0:
            A[i,i]=2+2*alpha+complex(0,dt*V[i]/hbar)
            B[i,i]=2-2*alpha-complex(0,dt*V[i]/hbar)
            A[i,i+1]=-alpha
            B[i,i+1]=alpha
            
    elif i==N-1:
            A[i,i]=2+2*alpha+complex(0,dt*V[i]/hbar)
            B[i,i]=2-2*alpha-complex(0,dt*V[i]/hbar)
            A[i,i-1]=-alpha
            B[i,i-1]=alpha
    else:
        A[i,i]=2+2*alpha+complex(0,dt*V[i]/hbar)
        B[i,i]=2-2*alpha-complex(0,dt*V[i]/hbar)
        A[i,i+1]=-alpha
        B[i,i+1]=alpha
        A[i,i-1]=-alpha
        B[i,i-1]=alpha
        
#print(A) 
#print("PSI",PSI)
#print(B)

########################################## Formulation of C
#print("PSI_2")
#print(PSI_2)
C=np.dot(np.linalg.inv(A),B)
# print("C")
# print(C)


# In[373]:


################################################ time forward

# for f in range(500):       #### for gif making purpose
for i in range(3*10**3):
    PSI=np.dot(C,PSI)

PSI_n2=(np.absolute(PSI))**2

PSI_n2=np.array(PSI_n2)*lencon*1e-6 #micrometer inverse
Xn=np.array(X)*1e6/lencon    #micrometer

##plotting

plt.plot(Xn,PSI_n2)
plt.ylabel("$||\psi(x)||^2$")
plt.xlabel("x ( $\mu$m )")
#plt.xlim(-25,25)

plt.axvline(x=-a/lencon, color="r")
plt.axvline(x=a/lencon, color="r")
#plt.savefig("snapshots1/snapshot_"+str(f)+".png")  ######for gif making purpose
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




