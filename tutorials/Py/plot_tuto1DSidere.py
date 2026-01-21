#
# Import section
#
import subprocess
import sys
import numpy as np
import fluidfoam
from pylab import figure, subplot, axis, xlabel, ylabel, show, savefig, plot
from pylab import title, matplotlib
import matplotlib.gridspec as gridspec
import matplotlib as mpl
#
# Change fontsize
#
matplotlib.rcParams.update({'font.size': 20})
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['lines.markersize'] = 10

#
# Change subplot sizes
#
gs = gridspec.GridSpec(2, 2)
gs.update(left=0.1, right=0.95, top=0.95,
          bottom=0.075, wspace=0.125, hspace=0.125)
#
# Figure size
#
figwidth = 8
figheight = 12

#########################################
#
# Loading experimental results
#
#
#
#
rho_s=2720.
rho_f=1000.

Htot=0.
#0.692
print ('rho_m^max=',0.39*(rho_s-rho_f)+rho_f)

#
# load exp
#
#rhom
npzfiles= np.load('DATA/Sidere_concentration_Sidc6.npz')
#
# DNS model 
#
rhom_sidc6_0d      = npzfiles['arr_0']
z_sidc6_0d         = npzfiles['arr_1']
rhom_sidc6_0_15d   = npzfiles['arr_2']
z_sidc6_0_15d      = npzfiles['arr_3']
rhom_sidc6_0_44d   = npzfiles['arr_4']
z_sidc6_0_44d      = npzfiles['arr_5']
rhom_sidc6_6d      = npzfiles['arr_6']
z_sidc6_6d         = npzfiles['arr_7']

#Ue
npzfiles= np.load('DATA/Sidere_pressure_Sidc6.npz')
#
# DNS model 
#
Ue_sidc6_0d      = npzfiles['arr_0']
z_Ue_sidc6_0d    = npzfiles['arr_1']
Ue_sidc6_0_15d   = npzfiles['arr_2']
z_Ue_sidc6_0_15d = npzfiles['arr_3']
Ue_sidc6_0_44d   = npzfiles['arr_4']
z_Ue_sidc6_0_44d = npzfiles['arr_5']
Ue_sidc6_6d      = npzfiles['arr_6']
z_Ue_sidc6_6d    = npzfiles['arr_7']

#########################################
# Loading OpenFoam results
#
#
#
#
case = '1DSidere'
basepath = '../laminar/'
sol = basepath + case + '/'

Nx = 1
Ny = 80
Nz = 1

eps_file = sol + case + '.eps'

#
# Reading SedFoam results
#
try:
    proc = subprocess.Popen(
        ['foamListTimes', '-case', sol], stdout=subprocess.PIPE)
except:
    print("foamListTimes : command not found")
    print("Do you have load OpenFoam environement?")
    sys.exit(0)
output = proc.stdout.read()
print(output)
#tread = list(output.split('\n'))
tread = output.decode().rstrip().split('\n')


del tread[-1]
Nt = len(tread)
print(Nt)
time = np.zeros(Nt)
X, Y, Z = fluidfoam.readmesh(sol)
alphat = np.zeros((Ny, Nt))
Ue = np.zeros((Ny, Nt))

k = -1
for t in tread:
    print(t)
    print("Reading time: %s s" % t)
    k = k + 1

    alphat[:, k] = fluidfoam.readscalar(sol, t + '/', 'alpha.a')
    Ue[:, k] = fluidfoam.readscalar(sol, t + '/', 'p_rbgh')
    time[k] = float(t)


#
# parameter
#
zmin = 0.
zmax = np.max(Y)

tmax = 518400.


#
# calcul zint et zint2
#
if Nt > 1:

    asint2 = 0.5*(0.39+0.28)
    asint  = 0.5*(0.39+0.)
    zint = np.zeros(Nt)
    zint2 = np.zeros(Nt)

    for i in np.arange(Nt):
        # zint
        toto = np.where(alphat[:, i] < asint)
        if np.size(toto) == 0:
            zint[i] = Y[Ny - 1]
        else:
            zint[i] = Y[toto[0][0]]
    # zint2
        toto2 = np.where(alphat[:, i] <= asint2)
        if np.size(toto2) == 0:
            zint2[i] = Y[0]
        else:
            zint2[i] = Y[toto2[0][0]]


figure(num=1, figsize=(figwidth*3/2, figheight),
           dpi=60, facecolor='w', edgecolor='w')
plot(time,zint,'-b',time,zint2,'-r')
ylabel('z (m)')
xlabel('t (s)')
axis([0 ,time[Nt-1], zmin, zmax])

#show()
#show(block=False)

#toto=raw_input('hit a key to continue')

# t=0 day
inum1=0

# t=0.15 day
toto=np.where(time[:]>=0.15*24.*3600.)
if (np.size(toto)==0):
    inum2=Nt-1
else:
    inum2=toto[0][0]
print (inum2,time[inum2],0.15*24.*3600.)

# t=0.44 day
toto=np.where(time[:]>=0.44*24.*3600.)
if (np.size(toto)==0):
    inum3=Nt-1
else:
    inum3=toto[0][0]
print (inum3,time[inum3],0.44*24.*3600.)

# t=6 day
toto=np.where(time[:]>=6.*24.*3600.)
if (np.size(toto)==0):
    inum4=Nt-1
else:
    inum4=toto[0][0]
print (inum4,time[inum4],6.*24.*3600.)


#
# Figure 2: Concentration profiles
#
figure(num=2, figsize=(figwidth, figheight),
           dpi=60, facecolor='w', edgecolor='w')

gs = gridspec.GridSpec(2, 1)
gs.update(left=0.075, right=0.965, top=0.94, bottom=0.12,wspace=0.3, hspace=0.25)

# Set the ticks and labels...
x_ticks = np.array([1200, 1400, 1600, 1800])
x_labels = ['1200', '1400', '1600', '1800']
y_ticks = np.array([-0.6, -0.4, -0.2, 0])
y_labels = ['-0.6', '-0.4', '-0.2', '0']

rhom_min=1450.
rhom_max=1800.
rho_leg=1680

ax13 = subplot(gs[0,0])

#
# 0 day
#
p1=ax13.plot(rhom_sidc6_0d,z_sidc6_0d,'--b',label='t=0 day')
p1n=ax13.plot(alphat[:,0]*(rho_s-rho_f)+rho_f,Y[:],'-r') 
ax13.annotate('t=0 day', xy=(alphat[int(Ny/2)-5,0]*(rho_s-rho_f)+rho_f, 
                                    Y[int(Ny/2)-5]),  xycoords='data',
                xytext=(rho_leg, .55), textcoords='data',
                arrowprops=dict(arrowstyle="->")
                )
#
# 0.15 day
#
p2=ax13.plot(rhom_sidc6_0_15d,z_sidc6_0_15d,'--b',label='t=0.15 day')
p2n=ax13.plot(alphat[:,inum2]*(rho_s-rho_f)+rho_f,Y,'-r') 
ax13.annotate('t=0.15 day', xy=(alphat[int(Ny/2)-17,inum2]*(rho_s-rho_f)+rho_f, 
                                       Y[int(Ny/2)-17]),  xycoords='data',
                xytext=(rho_leg, 0.475), textcoords='data',
                arrowprops=dict(arrowstyle="->")
                )
#
# 0.44 day
#
p3=ax13.plot(rhom_sidc6_0_44d,z_sidc6_0_44d,'--b',label='t=0.44 day')
p3n=ax13.plot(alphat[:,inum3]*(rho_s-rho_f)+rho_f,Y,'-r') 
ax13.annotate('t=0.44 day', xy=(alphat[int(Ny/2)-25,inum3]*(rho_s-rho_f)+rho_f, 
                                       Y[int(Ny/2)-25]),  xycoords='data',
                xytext=(rho_leg, 0.4), textcoords='data',
                arrowprops=dict(arrowstyle="->")
                )
#
# 6 days
#
p4=ax13.plot(rhom_sidc6_6d,z_sidc6_6d-Htot,'--b',label='t=6 days')
p4n=ax13.plot(alphat[:,inum4]*(rho_s-rho_f)+rho_f,Y,'-r') 
ax13.annotate('t=6 days', xy=(alphat[int(Ny/2)+8,inum4]*(rho_s-rho_f)+rho_f, 
                                     Y[int(Ny/2)+8]),  xycoords='data',
                xytext=(1550, Y[int(Ny/2)+23]), textcoords='data',
                arrowprops=dict(arrowstyle="->")
                )

title('SIDC 6')
axis([rhom_min, rhom_max, zmin, zmax])
ylabel('z (m)')
xlabel(r'$\rho_m$ (kg/m$^3$)')


Ue_max=3000.
# Set the ticks and labels...
x_ticks = np.array([0, 1000, 2000,3000])
x_labels = ['0', '1000','2000', '3000']

#
# Figure 4: Excess pore pressure profiles
#
ax23 = subplot(gs[1,0])

Ue_max=3000.
# Set the ticks and labels...
x_ticks = np.array([0, 1000, 2000,3000])
x_labels = ['0', '1000','2000', '3000']
#
# 0 day
#
p1=ax23.plot(Ue_sidc6_0d,z_Ue_sidc6_0d,'--b',label='t=0 day')
p1n=ax23.plot(Ue[:,1],Y[:],'-r') 
ax23.annotate('t=0 day', xy=(Ue[int(Ny/2)+10,1],Y[int(Ny/2)+10]),  xycoords='data',
                xytext=(2050., Y[int(Ny/2)+20]), textcoords='data',
                arrowprops=dict(arrowstyle="->")
                )
#
#
# 0.15 day
#
p2=ax23.plot(Ue_sidc6_0_15d,z_Ue_sidc6_0_15d,'--b',label='t=0.15 day')
p2n=ax23.plot(Ue[:,inum2],Y[:],'-r') 
ax23.annotate('t=0.15 day', xy=(Ue[int(Ny/2),inum2],Y[int(Ny/2)]),  xycoords='data',
                xytext=(2050., Y[int(Ny/2)+10]), textcoords='data',
                arrowprops=dict(arrowstyle="->")
                )
#
# 0.44 day
#
p3=ax23.plot(Ue_sidc6_0_44d,z_Ue_sidc6_0_44d,'--b',label='t=0.44 day')
p3n=ax23.plot(Ue[:,inum3],Y[:],'-r') 
ax23.annotate('t=0.44 day', xy=(Ue[int(Ny/2)-10,inum3],Y[int(Ny/2)-10]),  xycoords='data',
                xytext=(2050., Y[int(Ny/2)]), textcoords='data',
                arrowprops=dict(arrowstyle="->")
                )
#
#
# 6 days
#
p4=ax23.plot(Ue_sidc6_6d,z_Ue_sidc6_6d,'--b',label='t=6 days')
p4n=ax23.plot(Ue[:,inum4],Y[:],'-r') 
ax23.annotate('t=6 days', xy=(Ue[int(Ny/2)-35,inum4],Y[int(Ny/2)-35]),  xycoords='data',
                xytext=(2050., Y[int(Ny/2)-15]), textcoords='data',
                arrowprops=dict(arrowstyle="->")
                )

axis([0, Ue_max, zmin, zmax])
xlabel('Ue (Pa)')
ylabel('z (m)')

show()

savefig('Figures/res2_tuto5.png', facecolor='w', edgecolor='w', format='png')

show(block=False)
# Fix Python 2.x.
#try: input = raw_input
#except NameError: pass
#tto = input("Hit a key to close the figure")
