#Example of Abidoye et al (2014) case.
#Viscosity ratio = 500, Applied pressure drop = 15 kPa.

#Part - 1 (Code set-up)
#Import all essential PYTHON libraries
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#Part - 2 (input essential parameters)
#Input/ tune as per requirement the user defined parameters (fluid properties, porous medium properties)
theta = 40                           # contact angle
rp = 30e-6                        # pore body radii (average at macroscale)
rt = 12e-6                         # pore throat radii (average at macroscale)
beta = math.atan(-(rt/rp)+1)        # orientation angle
Rt_min = 8e-6                     # minimum pore throat radius @pore-scale
phi = 0.38                          # porosity
ift = 0.021                        # interfacial tension
t = 1000                            # number of saturation data points for analysis
swr = 0.138857                       # residual water saturation
mu_w = 8.9e-4                      # viscosity of wetting phase
mu_n = 0.445                       # viscoisty of non-wetting phase
delta_p = 10                        # applied pressure drop between inlet and outlet boundaries
L = 0.04                           # length between inlet and outlet boundaries where pressure drop is acting

#Part - 3 (Function definition to compute the grain radius, rg and the rugosity parameter, b) 
from scipy.optimize import minimize
def f(x):
    rg  = x[0]
    b = x[1]
    theta_radians  = math.radians(theta) 
    increament  = (1-swr)/t
    sw = np.zeros(t)
    Rt = np.zeros(t)
    Rt_m = np.zeros(t)
    sw[0] = swr
    for i in range(1,t):
        sw[i] = sw[i-1] + increament
    sn = 1-sw
    A = -((phi*sn*((1-math.tan(beta))**3)*(rg)**3)**(1/3))/3
    B1 = 8+(27*(1-math.tan(beta))**3)
    B2 = sn*(8+(27*(1-math.tan(beta))**3))
    B3 = phi*sn*(8+(27*(1-math.tan(beta))**3))
    B = B1/((B2-B3)**(4/3))
    DRtDsn = A*B
    Rt[0] = Rt_min
    Rt_m[0] = Rt_min
    for i in range(1,t):
        Rt[i] = Rt[i-1] - DRtDsn[i-1]*(increament)
        Rt_m[i] = Rt[i]/(((sw[i]**b)+sn[i]**2)*math.exp(sw[i]))
    return (((Rt_m[500]-rt)**2)/(rt))

#Part - 4 (Minimisation algorithm to compute the grain radius and the heuristic parameter)
res = minimize(f,[400e-6, 1.2],method='Nelder-Mead',tol=1e-30)
print(res)                          # Print the output of the results
rg = res.x[0]                       # Copy the value of the 'res' output from above to this variable 
b = res.x[1]

#Part - 5 (Compute new set of pore throat radii from new values of rg and a)
theta_radians  = math.radians(theta)
increament  = (1-swr)/t
sw = np.zeros(t)
Rt = np.zeros(t)
Rt_m = np.zeros(t)
sw[0] = swr
for i in range(1,t):
    sw[i] = sw[i-1] + increament
sn = 1-sw
A = -((phi*sn*((1-math.tan(beta))**3)*(rg)**3)**(1/3))/3
B1 = 8+(27*(1-math.tan(beta))**3)
B2 = sn*(8+(27*(1-math.tan(beta))**3))
B3 = phi*sn*(8+(27*(1-math.tan(beta))**3))
B = B1/((B2-B3)**(4/3))
DRtDsn = A*B
Rt[0] = Rt_min
Rt_m[0] = Rt_min
for i in range(1,t):
    Rt[i] = Rt[i-1] - DRtDsn[i-1]*(increament)
    Rt_m[i] = Rt[i]/(((sw[i]**b)+sn[i]**2)*math.exp(sw[i]))

#Part - 6 (Compute the static capillary pressure and compare with literature data)
pc = (2*(ift*math.cos(theta_radians)))/(Rt_m*1000)          # in kPa
pc1 = pd.read_excel(r'C:\Users\saideep.pavuluri\Desktop\dynPcSdp.xlsx',sheet_name='Das_M500', skiprows=1, usecols="A,A:H")
#pc1.head()
plt.plot(sw,pc,color='b',linewidth=3,label='Semi-analytical model')
plt.scatter(pc1['sw'],pc1['pc_static'],color='k',label='Experiment')
plt.xlim(0, 1)
plt.ylim(0, 5)
plt.legend()
plt.savefig('pcs_M500_p10.pdf',format='pdf')
plt.show()  # Displays the static capillary pressure

#Part - 7 (Compute the change in permeability as a function of saturation and 'DRtDsw')
Asw = ((phi*sw*((1-math.tan(beta))**3)*(rg)**3)**(1/3))/3
B1sw = 8+(27*(1-math.tan(beta))**3)
B2sw = sw*(8+(27*(1-math.tan(beta))**3))
B3sw = phi*sw*(8+(27*(1-math.tan(beta))**3))
Bsw = B1sw/((B2sw-B3sw)**(4/3))
DRtDsw = Asw*Bsw
DkDsw = (DRtDsw*Rt_m)*phi   # Eq. 9 in manuscript
DkDsn = (-1*DRtDsn*Rt_m)*phi

#Part - 8 (compute the Darcy flow velocity)
u = ((delta_p-pc)/(pc))*(ift*math.cos(theta_radians)*rt*phi)/(4*L*(mu_n*sn+mu_w*sw)) # Eq. 14 in manuscript

#Part - 9 (compute viscous ressure drops of wetting and non-wetting fluids)
delta_pw = u*mu_w*L*(1/DkDsw)   # Eq. 15 in manuscript
delta_pn = u*mu_n*L*(1/DkDsn)   # Eq. 15 in manuscript

#Part - 10 (compute the dynamic capillary pressure)
pc_dynamic = pc + (delta_pn - delta_pw)/1000    # Eq. 18 in manuscript

#Part - 11 (Plot the dynamic capillary pressure and compare with literature data)
plt.plot(sw,pc_dynamic,color='r',linewidth=3,label='Semi-analytical model')
plt.scatter(pc1['sw_10'],pc1['pc_10'],color='k',label='Experiment')
plt.xlim(0.2, 1)
plt.ylim(-2, 10)
plt.legend()
#plt.savefig('pcd_M500_p10.pdf',format='pdf')   # Uncomment to save plot
plt.show()  # Display the result
