#A sample model in which essential data has to be populated by the user.

#Part - 1 (Code set-up)
#Import all essential PYTHON libraries
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#Part - 2 (input essential parameters)
#Input/ tune as per requirement the user defined parameters (fluid properties, porous medium properties)
theta =                            # contact angle
rp =                         # pore body radii (average at macroscale)
rt =                          # pore throat radii (average at macroscale)
beta = math.atan(-(rt/rp)+1)        # orientation angle
Rt_min =                      # minimum pore throat radius @pore-scale
phi =                           # porosity
ift =                         # interfacial tension
t =                             # number of saturation data points for analysis
swr =                        # residual water saturation
mu_w =                        # viscosity of wetting phase
mu_n =                        # viscoisty of non-wetting phase
delta_p =                         # applied pressure drop between inlet and outlet boundaries
L =                            # length between inlet and outlet boundaries where pressure drop is acting

#Part - 3 (Function definition to compute the grain radius, rg and the rugosity parameter, b) 
from scipy.optimize import minimize
def f(x):
    rg  = x[0]                      # grain radius
    a = x[1]                        # rugosity (heuristic) parameter
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
    DRtDsn = A*B                    # Refer to Eq. 3 in the manuscript
    Rt[0] = Rt_min
    Rt_m[0] = Rt_min                
    for i in range(1,t):
        Rt[i] = Rt[i-1] - DRtDsn[i-1]*(increament)  # NOTE: (-) for variable 'A' above and (-) for DRtDsn here represent Eq. 2 in manuscript 
        Rt_m[i] = Rt[i]/(((sw[i]**a)+sn[i]**2)*math.exp(sw[i])) # Eq. 4 in manuscript
    return (((Rt_m[500]-rt)**2)/(rt))                           # return the function output

#Part - 4 (Minimisation algorithm to compute the grain radius and the heuristic parameter)
res = minimize(f,[400e-6, 1.2],method='Nelder-Mead',tol=1e-30)
#print(res)                          # Print the output of the results
rg = res.x[0]                       # Copy the value of the 'res' output from above to this variable 
a = res.x[1]

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
    Rt_m[i] = Rt[i]/(((sw[i]**a)+sn[i]**2)*math.exp(sw[i]))

#Part - 6 (Compute the static capillary pressure and compare with literature data)
pc = (2*(ift*math.cos(theta_radians)))/(Rt_m*1000)          # in kPa
pc1 = pd.read_excel(r'C:\Users\saideep.pavuluri\Desktop\dynPcSdp.xlsx',sheet_name='Vahid_M1', skiprows=1, usecols="A,A:H") #literature data. Modify as per your file location
#pc1.head()
plt.plot(sw,pc,color='b',linewidth=3,label='Semi-analytical model')
plt.scatter(pc1['sw'],pc1['pc_static'],color='k',label='Static PNM')
plt.xlim(0, 1)
#plt.ylim(0, 40)
plt.legend()
plt.show()                      # Displays the static capillary pressure

#Part - 7 (Compute the change in permeability as a function of saturation and 'DRtDsw')
Asw = ((phi*sw*((1-math.tan(beta))**3)*(rg)**3)**(1/3))/3
B1sw = 8+(27*(1-math.tan(beta))**3)
B2sw = sw*(8+(27*(1-math.tan(beta))**3))
B3sw = phi*sw*(8+(27*(1-math.tan(beta))**3))
Bsw = B1sw/((B2sw-B3sw)**(4/3))
DRtDsw = Asw*Bsw
k =                             #Variable describing heterogeniety effects and impact of rel perm for the wetting phase
j =                             #Variable describing heterogeniety effects and impact of rel perm for the non-wetting phase
DkDsw = (2*DRtDsw*Rt_m)*phi*k   # Eq. 9 in manuscript
DkDsn = -(2*DRtDsn*Rt_m)*phi*j

#Part - 8 (compute the Darcy flow velocity)
u = ((delta_p-pc)/(pc))*(ift*math.cos(theta_radians)*rt*phi)/(4*L*(mu_n*sn+mu_w*sw)) # Eq. 14 in manuscript

#Part - 9 (compute viscous ressure drops of wetting and non-wetting fluids)
delta_pw = u*mu_w*L*(1/DkDsw)   # Eq. 15 in manuscript
delta_pn = u*mu_n*L*(1/DkDsn)   # Eq. 15 in manuscript

#Part - 10 (compute the dynamic capillary pressure)
pc_dynamic = pc + (delta_pn + delta_pw)/1000    # Eq. 17 in manuscript

#Part - 11 (Plot the dynamic capillary pressure and compare with literature data)
plt.plot(sw,pc_dynamic,color='r',linewidth=3,label='Semi-analytical model')
plt.scatter(pc1['sw'],pc1['pc'],color='k',label='Dynamic PNM')
plt.xlim(0.2, 1)
#plt.ylim(0, 40)
plt.legend()
#plt.savefig('name.pdf',format='pdf')     # Uncomment to save plot
plt.show()                      # Display the result
