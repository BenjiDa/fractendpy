import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sga import shearonplane, pole, principalstress, stcoordline
import pdb

trend_s1 = 26 #Stress orientations from Boyle and Zoback 2014
plunge_s1 = 0
trend_s3 = 116

ncontours = 20

sigma1 = 50       
sigma2 = 20     
sigma3 = 5

sorted_sigma = [sigma1, sigma2, sigma3]
sigmad = sigma1 - sigma3
stress_tensor = np.diag(sorted_sigma)


# convert stress orientations into radians
trend_s1_rad = np.radians(trend_s1)# * np.pi / 180
plunge_s1_rad = np.radians(plunge_s1)# * np.pi / 180 
trend_s3_rad = np.radians(trend_s3)# * np.pi / 180


#   define a 3d space
# for all directions in 3-space, calculate normal and shear stresses 
#   on all 3d surfaces    
phi_index = 0  
theta_index = 0  
increment = 10
phi_min = 90  
phi_max = 180  
theta_min = 0  
theta_max = 360  
phi_n = (phi_max - phi_min)//increment + 1  #divide to make int
theta_n = (theta_max - theta_min)//increment + 1 
sigmaN = np.zeros((theta_n, phi_n))  
tau = np.zeros((theta_n, phi_n)) 


for idp, phi in enumerate(range(0,91,increment)):

    phi_rad = np.radians(phi)# * np.pi/180

    for idt, theta in enumerate(range(0,361,increment)):

        theta_rad = np.radians(theta)# * np.pi / 180

        #   convert pole to strike and dip 
        [strike, dip] = pole(theta_rad, phi_rad, 0) 

        #   calculate normal and shear stress on the plane
        [stress_fracture, dc_stress_fracture, R] = shearonplane(stress_tensor, trend_s1_rad, plunge_s1_rad, trend_s3_rad, strike, dip) 
        print('Dip: ', dip, 'Stress :', stress_fracture)

        #np.nan_to_num(stress_fracture, copy=False, nan=0.0) # This line converts nan values into zero, nan values are generated when
        #the dip of the surface is parallel to sigma 1 and sigma 3

        #   save normal and shear stresses for later calculation
        sigmaN[idt][idp] = stress_fracture[0][0]#(1,1) 
        tau[idt][idp] = stress_fracture[2][0]#(3,1)



deltaP = increment * np.pi / 180 
phiP = np.arange(np.pi/2, np.pi+deltaP, deltaP)
phiP = phiP - np.pi/2  
thetaP = np.arange(0, 2*np.pi+deltaP, deltaP)  
[phiP, thetaP] = np.meshgrid(phiP, thetaP)

#   equal area projection 
dp = np.sqrt(1 - np.sin(phiP)) 
xeqarea = dp * np.sin(thetaP)  
yeqarea = dp * np.cos(thetaP)   
rPrim = 1  
xPrim = np.arange(rPrim*-1, rPrim, 0.0001) 
yPrim = np.sqrt(np.square(rPrim) - np.square(xPrim)) 



#   convert principal stress orientations into cartesian coords for equal
#   area
[pstress, dCp] = principalstress(stress_tensor, trend_s1_rad, plunge_s1_rad, trend_s3_rad) 
trendS2rad = pstress[1,1] 
plungeS2rad = pstress[1,2] 
plungeS3rad = pstress[2,2] 
[xS1, yS1] = stcoordline(trend_s1_rad, plunge_s1_rad, 1)   
[xS2, yS2] = stcoordline(trendS2rad, plungeS2rad, 1)   
[xS3, yS3] = stcoordline(trend_s3_rad, plungeS3rad, 1)



fig, ax = plt.subplots(1, figsize=(6,6))

sn = plt.contourf(xeqarea, yeqarea, tau, ncontours)#, 'EdgeColor', 'none')  

ax.plot(xPrim, yPrim, '-k', linewidth=1)  
ax.plot(xPrim, yPrim*-1, '-k', linewidth=1) 
#ax.plot(xFractures, yFractures, '.r',markeredgecolor='k', markersize=15) 
ax.plot(xS1, yS1, 's', markersize=10, markeredgecolor='k', markerfacecolor='w')
ax.annotate(text='\u03C3'+'1', xy=(xS1+0.05, yS1+0.05))  
ax.plot(xS2, yS2, 'd', markersize=10, markeredgecolor='k', markerfacecolor='w') 
ax.annotate(text='\u03C3'+'2', xy=(xS2+0.05, yS2+0.05))  
ax.plot(xS3, yS3, '^', markersize=10, markeredgecolor='k', markerfacecolor='w')
ax.annotate(text='\u03C3'+'3', xy=(xS3+0.05, yS3+0.05)) 
ax.set_xticks([])
ax.set_yticks([])

plt.gca().set_aspect('equal', adjustable='box')

plt.xlim([-1.05, 1.05])  
plt.ylim([-1.05, 1.05])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

# creating colorbar
fig.colorbar(sn, cax = cax, orientation = "vertical")

#fig.axes[0].set_title('{}'.format(dataset_name) + ' n={}'.format(n_fractures))

fig.savefig("./{}_stereo_test.pdf".format('tau'), dpi=300, transparent=True)
plt.show()
