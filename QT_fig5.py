#!/usr/bin/env python
# coding: utf-8

# In[14]:


################# Reproducing Fig5 (varry width)  ##############


import numpy as np
from math import exp, pi, pow,sqrt
from matplotlib import pyplot as plt



##################################################self decleared values (IN SI UNITS)

x0=-14*1e-6 #meter
k0=2*1e6   #meter inverse
sigmak = 0.7*1e6;
#sigmak=0.5 #micrometer^-1
sigmax=1/(2*sigmak)
hbar = 1.054*1e-34
m = 9.1*1e-31
E0 = hbar**2*k0**2/(2*m)
dt = 0.05*1e-9 




######################################## conversion factors from si to natural

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

###########################################################   Wave function initialization
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
X=np.arange(-xwidth,xwidthR,xgrid) #natural unit
print('xwidth=',xwidth)
print('xgrid',xgrid)




# #print(PSI_2)
# plt.plot(X,PSI_2,label="PSI_2 in natural space")
# plt.legend()
# plt.show()


# PSI_c2=np.array(PSI_2)*lencon*1e-6 #micrometer inverse
# Xn=np.array(X)*1e6/lencon    #micrometer
# plt.plot(Xn,PSI_c2,label="PSI_2 in x space")
# plt.ylabel("$||\psi(x)||^2$")
# plt.xlabel("x($\mu$m)")
# plt.legend()
# plt.show()



##############################################       potential initialization
v0list = [1.5]

v0 = v0list[0]
a = 1e-6  #micrometer
alist =np.arange(0.0,2.000001,0.05)*1e-6*lencon
P=[]
for a in alist:
    PSI = [psi(x) for x in X]
    V=[]
    for x1 in X:
        if abs(x1)>a:
            v=0*econ
        else:
            #v=0
            v=v0*E0
        V.append(v)

####################################################### A and B and C matrix formulation
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

##################################################################### Formulation of C
    #print("PSI_2")
    #print(PSI_2)
    C=np.dot(np.linalg.inv(A),B)
    # print("C")
    # print(C)
###################################################################### time forward

 # for f in range(500):
    for i in range(30*10**2):
        PSI=np.dot(C,PSI)


    PSI_2V1=(np.absolute(PSI))**2
    PSI_2V1=np.array(PSI_2V1)*lencon*1e-4 #meter inv
    X_m=np.array(X)*1e5/lencon     #meter
    P.append(list(PSI_2V1).index(max(PSI_2V1[int(N/2):])))
#     plt.xlim(0,7)
#     plt.ylim(top=3.5)
#     plt.plot(X_m,PSI_2V1,label='trans wave for a='+str(2*a/lencon))
# plt.ylabel("$|\psi|^2 x10^4$")
# plt.xlabel("x (m) x$10^{-5}$")
# plt.legend()
# plt.show()
print(P)


# In[15]:


#######################       plotting

shiftX=[]
for p in P:
    shiftX.append((X_m[p]-X_m[P[0]])*10)
plt.plot(2*alist*1e6/lencon,shiftX,"-p")
plt.xlabel("Barrier width ($\mu$m)")
plt.ylabel("Shift of the peak ($\mu$m)")


# In[ ]:





# In[ ]:





# In[ ]:




