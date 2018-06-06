
# coding: utf-8

# In[464]:


import matplotlib.pyplot as plt
import numpy as np
import os
import math
import scipy.interpolate
from scipy import ndimage
import pylab as pl
import sys
import pandas as pd
import h5py

data = 'HM_vorticity' # [Potential, Vorticity, HM_vorticity] Choose the data that you want to visualize

current_dir2 = os.getcwd()
current_dir = 'U:\\fortran\\code_HM'
result_Field= current_dir + "\%s.dat" %data
print("Current using data is = ", result_Field)
Field = np.fromfile(result_Field, dtype=float)
print("Number of elements in Field.dat =",len(Field))


# In[465]:


print('Python version ' + sys.version)
print('Pandas version ' + pd.__version__)

So, The time is given by [dt, nstep, noutt, noutf]
And number of elements in field-2D data is given by value of [nx, ny] in mod_global.f03
When nx and ny is given as power(n) of 2, then number of field is (nx+ny) power of 2
# ### Simulation parameters

# In[466]:


# Grid
nx = 256
ny = 256
Nx = int(math.log2(nx))
Ny = Nx
Field_cycle = 2**(Nx+Ny)
deltax = 2 * np.pi / nx
# Time
dt = 0.0001
nstep = 3000
noutt = 10
noutf = 100
Nb_time = int( (nstep / noutf) + 1 ) 
t_tot = dt * nstep

# Other parameters
omega = 0.66
visc = 3e-6

print('Number of grid (nx, ny) =', '(',nx,ny,')')
print('Total time t =',t_tot,'(s)')
print( 'delta t = ', dt , '/', '(delta x)**2 =', (deltax)**2)
print('Filtering omega =',omega,'/','Viscosity =',visc)
print('delta x = ', deltax)


# In[467]:


time = []
for p in range(Nb_time):
    #print(Field[p*(Field_cycle+1)])
    time.append(Field[p*(Field_cycle+1)])
#print(time)


# In[468]:


Field1=[[[] for i in range(nx)] for j in range(Nb_time)]
for p in range(Nb_time):
    for q in range(ny):
        for r in range(nx):
            Field1[p][q].append(Field[p*(Field_cycle+1) + q * ny + (r + 1)])
            
# In this Field table, we have a (Nb_time) dimensions, (nx) dimensions, (ny) dimensions.
# So in each iterations, we fix the value of time and then value of nx and then inside we save all the value of ny when
# nx is fixed. We repeat this cycle by changing value of nx and value of time.


# In[469]:


'''
for p in range(Nb_time):
    print("time is, t=",p)
    for q in range(nx):
        print("Ny is =", q)
        print("Field1(Nx,Ny)")
        for r in range(nx):
            print("Field1(",r,q,") =",Field1[p][q][r])
'''


# ### Contour graph for field data

# In[470]:


file_path = "U:\\fortran\\code_HM\\results\\%s\\nx%sny%st%somega%svisc%s" %(data,nx,ny,t_tot,omega,visc)
if not os.path.exists(file_path):
    os.makedirs(file_path)   


# In[471]:


for i in range(Nb_time):
    fig, axes = pl.subplots(dpi=100)
    c1 = axes.contourf(Field1[i])
    pl.colorbar(c1, ax=axes);
    # pl.colorbar(c1, ax=axes, cmap='jet');
    # axes.set_title(r'$\nabla^2 \phi - \phi$, time = %i' %i)
    axes.set_title('%s, time = %i' %(data,i))
    fig.savefig(file_path + "\%s_t%i.png" %(data,i))


# In[472]:


# Interpolated images
xi = np.linspace(0, 2*np.pi, nx)
yi = np.linspace(0, 2*np.pi, ny)
xi, yi = np.meshgrid(xi, yi)


# ### This part is for creating 'Time evolution' animation

# In[473]:


import re

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split('(\d+)', text) ]


# In[474]:


import os
filepath = 'U:\\fortran\\code_HM\\results\\%s\\nx%sny%st%somega%svisc%s' %(data,nx,ny,t_tot,omega,visc)
filenames = os.listdir('U:\\fortran\\code_HM\\results\\%s\\nx%sny%st%somega%svisc%s' %(data,nx,ny,t_tot,omega,visc))
filenames.sort(key=natural_keys)
# print(filenames)
# We need to make them in increasing order in order to create the animation


# In[475]:


import imageio
images = []
for filename in filenames:
    images.append(imageio.imread('U:\\fortran\\code_HM\\results\\%s\\nx%sny%st%somega%svisc%s\\' %(data,nx,ny,t_tot,omega,visc) +filename ))
imageio.mimsave('U:\\fortran\\code_HM\\results\\%s\\nx%sny%st%somega%svisc%s\\time_evolution.gif' %(data,nx,ny,t_tot,omega,visc), images)

