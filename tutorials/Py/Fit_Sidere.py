import numpy as np
from pylab import *

from scipy import optimize
from scipy import signal
#
# Change fontsize
#
matplotlib.rcParams.update({'font.size': 20})

#
#
#
lw=2.
ms=8.

#
# Effective stress
#
fileName='DATA/Sidere_EffectiveStress.npz'

npzfiles= np.load(fileName)

Sigma_sidc2 = npzfiles['arr_0']
e_sidc2     = npzfiles['arr_1']
Sigma_sidc3 = npzfiles['arr_2']
e_sidc3     = npzfiles['arr_3']
Sigma_sidc5 = npzfiles['arr_4']
e_sidc5     = npzfiles['arr_5']
Sigma_sidc6 = npzfiles['arr_6']
e_sidc6     = npzfiles['arr_7']

#
# Plot
#
Nc=500
c=linspace(0.01,0.6,Nc)

#%
#% Constante pour le regime de permeabilite
#%
# obtenu par best-fit sur les donnees de permeabilite
# sidc2-3 
B = - 7.93
# sidc5-6
#B = - 14.35

n = 3.+2./B
print ('dimension fractale: n=',n)
#n=2.85
#B=2./(n-3)

#%
#% Vitesse de chute d'un floc seul
#%
# sidc2-3 
Ws0=1.15e-4
# sidc5-6 
#Ws0=1e-3

#
# constante de la loi de contrainte effective
#
# sidc2-3
sigma0=75
sigma0_1=1.5e6
# sidc5-6 
#sigma0=5.
#sigma0_1=3.e8



rho_p=2720.
rho_f=1000.

#%
#% Concentration de gel
#%
c_gel=0.28
Cgel=c_gel*rho_p

#
# Effective stress closure
#
#cmax=0.39
cmax=0.45

sigma=sigma0*((1.-(c-c_gel)/cmax)**(B)-1.)
sigma1=sigma0_1*c**(-B)
sigma2=1.6e10*c**(19)

#
# 
#
rs=10e-6
xi=4.5
E=3e9
nu=0e0
mue=E/(2e0*(1e0+nu))
m=2e0/(9e0*sqrt(3e0))*mue*(2e0*rs)**2/(1e0-nu)

cmax2=0.9
sigma3=m/(pi*(2*rs)**2)*c*(c-c_gel)**(4.5)*3e0*(1e0+sin(0.5e0*pi*(2e0*(c-c_gel)/(cmax2-c_gel)-1e0)))

#clf()
figure(num=1, figsize=(18, 8),dpi=60, facecolor='w', edgecolor='w')

x_ticks = np.array([0.2, 0.3, 0.4])
x_labels = ['0.2', '0.3', '0.4']
y_ticks = np.array([0, 1000, 2000, 3000])
y_labels = ['0', '1e3','2e3', '3e3']


alphas_sidc2=1.-e_sidc2/(e_sidc2+1.)
alphas_sidc3=1.-e_sidc3/(e_sidc3+1.)
alphas_sidc5=1.-e_sidc5/(e_sidc5+1.)
alphas_sidc6=1.-e_sidc6/(e_sidc6+1.)

ax1=subplot(1,2,1)
title('Effective Stress')
p11 = ax1.plot(alphas_sidc2,Sigma_sidc2*1e3,'+k',label="sidc2",markersize=ms)
p12 = ax1.plot(alphas_sidc3,Sigma_sidc3*1e3,'xk',label="sidc3",markersize=ms)
p13 = ax1.plot(alphas_sidc5,Sigma_sidc5*1e3,'*k',label="sidc5",markersize=ms)
p14 = ax1.plot(alphas_sidc6,Sigma_sidc6*1e3,'sk',label="sidc6",markersize=ms)

p15=ax1.plot(c,sigma,label='Closure law',linewidth=lw)

#p16= ax1.plot(c,sigma2,'--m',label='Winterwerp',linewidth=lw)

#p16=ax1.plot(c,sigma1,'--b',label='Closure law bis',linewidth=lw)
#p18=ax1.plot(c,sigma3,label='Closure law 3',linewidth=lw)

handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles, labels,loc=2,prop={'size':16})
xlabel(r"$\alpha$")
ylabel(r"$\sigma_e$ (Pa)")
ax1.axis([0.2,0.45,0.,3e3])
xticks(x_ticks, x_labels)
yticks(y_ticks, y_labels)


#
# Permeability
#
fileName='DATA/Sidere_Permeability.npz'

npzfiles= np.load(fileName)

K_sidc2 = npzfiles['arr_0']
e_sidc2     = npzfiles['arr_1']
K_sidc3 = npzfiles['arr_2']
e_sidc3     = npzfiles['arr_3']
K_sidc5 = npzfiles['arr_4']
e_sidc5     = npzfiles['arr_5']
K_sidc6 = npzfiles['arr_6']
e_sidc6     = npzfiles['arr_7']

# Plot
alphas_sidc2=1.-e_sidc2/(e_sidc2+1.)
alphas_sidc3=1.-e_sidc3/(e_sidc3+1.)
alphas_sidc5=1.-e_sidc5/(e_sidc5+1.)
alphas_sidc6=1.-e_sidc6/(e_sidc6+1.)

#
# Loi de permeabilite
#
#
# Determination de B et n a partir d'un fit des donnees expe
#
alphas=np.concatenate((alphas_sidc2,alphas_sidc3,alphas_sidc5,alphas_sidc6),axis=0)
K=np.concatenate((K_sidc2,K_sidc3,K_sidc5,K_sidc6),axis=0)

#toto=np.where(alphas>=0.35)
#toto=np.where(alphas>0.28)
toto=np.where(alphas<0.35)

alphas_keep=alphas[toto]
K_keep=K[toto]

logx=log(alphas_keep)
y=log(K_keep)

# define our (line) fitting function
fitfunc = lambda p, logx: p[0] + p[1] * logx
errfunc = lambda p, logx, y: y - fitfunc(p, logx)
# 
pinit = [0.01, 0.01]
pfinal, success = optimize.leastsq(errfunc, pinit[:], args=(logx, y))
#
print ('B=',pfinal[1],' et n=',3.+2./B)

#%
#% Constante pour le regime entrave
#%
#% rho_floc impose => phi_max calcule
#%rho_floc = ...;
#%phi_max = c_gel*(rho_p-rho_f)/(rho_floc-rho_f);

#% phi_max impose => rho_floc calcule
phi_max = 0.85
rho_floc = rho_f + (rho_p-rho_f) * c_gel/phi_max
phi = (rho_p-rho_f)/(rho_floc-rho_f)*c

print ('phi_max=',phi_max)

Wsh= Ws0*(1.-c[c<=c_gel])**(n/2.)*(1.-phi[c<=c_gel])**(n/2.-1.)*(1.-phi[c<=c_gel]/phi_max)**phi_max
ch=c[c<=c_gel]

kh=(1.-ch)*rho_f/(ch*(rho_p-rho_f))*Wsh;

#
test=0
for i in range(np.size(ch)-1):
    dlnwshdlnc=(log(Wsh[i])-log(Wsh[i-1]))/(log(c[i])-log(c[i-1]));
    if ((dlnwshdlnc<=B+1) and (test==0)):
        xsi=c_gel/(0.5*(c[i]+c[i-1]));
        A=0.5*(Wsh[i]+Wsh[i-1])/(rho_p/rho_f-1.)*(xsi/c_gel)**(B+1);
        test=1;


Ws=np.zeros(Nc)

for i in range(np.size(c)-1):
    if (   (c[i]>=(c_gel/xsi) ) and  ( c[i-1]<=(c_gel/xsi) )   ):
            Ws_gel = Wsh[i]
            Ws[i] = Ws_gel
    elif (c[i]>c_gel/xsi):
        Ws[i] = Ws_gel*(xsi*c[i]/c_gel)**(B+1.)
    else:
        Ws[i] = Wsh[i];

k=(1.-c)*rho_f/(c*(rho_p-rho_f))*Ws

print ('Ws0=',Ws0)
print ('Ws_gel=',Ws_gel)
print ('xsi=',xsi)
print ('dimfract=',n)
print ('rho_floc=',rho_floc)
print ('sigma0=',sigma0)

Ws = Ws_gel*(xsi*c/c_gel)**(B+1.)
kp=(1.-c)*rho_f/(c*(rho_p-rho_f))*Ws

kx3=k*3
#kwinterwerp=1e-15*c**(-19)
kwinterwerp=1.6e-10*c**(-19)
#
# plot
#
y_ticks = np.array([0, 2e-6, 4e-6, 6e-6, 8e-6, 10e-6])
y_labels = ['0', '2e-6', '4e-6', '6e-6', '8e-6', '10e-6']

ax4=subplot(1,2,2)
title('Permeability')

p45= ax4.plot(ch,kh,'-r',label='Hindered settling',linewidth=lw)
p46= ax4.plot(c,kp,'-.g',label='Permeability',linewidth=lw)
p47= ax4.plot(c,k,'--b',label='K',linewidth=lw)

#p48= ax4.plot(c,kx3,'-b',label='Kx3',linewidth=lw)

p41 = ax4.plot(alphas_sidc2,K_sidc2,'+k',label="sidc2",markersize=ms)
p42 = ax4.plot(alphas_sidc3,K_sidc3,'xk',label="sidc3",markersize=ms)
p43 = ax4.plot(alphas_sidc5,K_sidc5,'*k',label="sidc5",markersize=ms)
p44 = ax4.plot(alphas_sidc6,K_sidc6,'sk',label="sidc6",markersize=ms)

#p49= ax4.plot(c,kwinterwerp,'--m',label='Winterwerp',linewidth=lw)

handles, labels = ax4.get_legend_handles_labels()
ax4.legend(handles[0:3], labels[0:3],loc=1,prop={'size':16})
ax4.set_yscale('log')
ax4.set_xscale('log')
xlabel(r'$\alpha$')
ylabel("K (m/s)")
ax4.axis([1e-1,1e0,1e-8,1e-2])


#savefig('../Results/FitSidereTHESISpres.eps', facecolor='w', edgecolor='w', format='eps')

show()        
