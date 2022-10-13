# /*%   FracTend.m - script to plot slip and dilatation tendency  
# %   
# %   Equations & code from:
# %       Morris et al., 1996 Geology
# %       Ferrill et al., 1999 GSA Today 
# %       Streit & Hillis, 2004 Energy
# %       Jolly & Sanderson, 1997 Journal of Structural Geology 
# %       Allmendinger et al., 2012 Structural Geology Algorithms, Cambridge
# %       University Press
# %
# %   David Healy & Tara Stephens  
# %   July 2018 
# %   d.healy@abdn.ac.uk*/
'''
    fractendpy converted to python by Ben Melosh

    fractendpy requires poles, 
'''

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from utils import stcoordline
from geometric import pole
from stress import principalstress
from shear import shearonplane

#from plot_mohr_oa import plot_mohr_oa

import pdb

## User inputs:

#   read in poles to specific fractures; tab-delimited text file, 
#    formatted as plunge then trend  
fractures_file = 'Utah_OA-Sills.txt' ;  # 'Utah_Thrusts_TS-RW. Utah_Sills_ArcRW. Utah_DefBands2RW. Utah_OA-Sills

fracture_poles = np.loadtxt(fractures_file)



#   read in stress magnitudes 
#   principal stresses in MPa 
sigma1 = 43       
sigma2 = 38.86     
sigma3 = 25 
# sigma1 = 100
# sigma2 = 48
# sigma3 = 29

Pf = 37.24 # pore fluid pressure in MPa, for fracture Opening Angle calculations


#   read in stress orientation 
#   e.g. for old case of SHmax azimuth of 135 
# normal fault system: Trend s1=0, Plunge s1=90 - use s3-trend to change orientation of stress field with s1 as the rotation axis
# thrust fault: Trend of s3 must be >90 from s1; Ps1=0
# strike-slip: Trend of s3 must be 90 from Trend of s1
trend_s1 = 68 
plunge_s1 =  3    
trend_s3 = 265 
# trend_s1 = 280 
# plunge_s1 = 60
# trend_s3 = 20


#   coefficient of friction & cohesion 
mu_static = 0.6 
cohesion = 0
sigmaN_mohr = 100
tau_mohr = cohesion + mu_static * sigmaN_mohr 

### End user inputs



# Convert fracture poles to radians and count number of fractures
n_fractures = len(fracture_poles)
fracture_poles_rad = fracture_poles * np.pi / 180  

sorted_sigma = [sigma1, sigma2, sigma3] 
sigmad = sigma1 - sigma3
stress_tensor = np.diag(sorted_sigma) #   Generate stress tensor. note the implicit convention: x//s1, y//s2, z//s3


    

# convert stress orientations into radians
trend_s1_rad = trend_s1 * np.pi / 180
plunge_s1_rad = plunge_s1 * np.pi / 180 
trend_s3_rad = trend_s3 * np.pi / 180




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


for idp, phi in enumerate(range(90,181,increment)):

    phi_rad = (phi - 90) * np.pi/180

    for idt, theta in enumerate(range(0,361,increment)):

        theta_rad = theta * np.pi / 180

        #   convert pole to strike and dip 
        [strike, dip] = pole(theta_rad, phi_rad, 0) 

        #   calculate normal and shear stress on the plane
        [stress_fracture, dc_stress_fracture, R] = shearonplane(stress_tensor, trend_s1_rad, plunge_s1_rad, trend_s3_rad, strike, dip) 

        #   save normal and shear stresses for later calculation
        sigmaN[idt][idp] = stress_fracture[0][0]#(1,1) 
        tau[idt][idp] = stress_fracture[2][0]#(3,1)

        


#   calculate tendencies - slip, dilatation and frac. suscep
#   calculate normalised slip tendency (Morris et al., 1996)

TsMax = (tau / sigmaN ).max() 
Ts = ( tau / sigmaN ) / TsMax  

#   calculate dilatation tendency (Ferril et al., 1999)
Td =  ( sorted_sigma[0] - sigmaN ) / ( sorted_sigma[0] - sorted_sigma[2] )   

#   calculate shear stress/ dilation displacement ratio 
#   from Delaney et al.(1988)
TD = tau / (sorted_sigma[2] + Pf) 

#   calculate fracture susceptibility
Sf = sigmaN - ( tau / mu_static )     

#   calculate muOA (opening angle), Jolly & Sanderson, 1997
OA = tau / (Pf - sigmaN)
muOAfracture = np.arctan(OA) * (180/np.pi)




#   Calculate and display various stress ratios:
Phi = ( sigma2 - sigma3 ) / ( sigma1 - sigma3 )  
R = (sigma1 - sigma2) / (sigma1 - sigma3) 
Rprime = ( Pf - sigma3 ) / ( sigma1 - sigma3 )  




#   plot azimuthal variation of tendencies, with poles to fractures
#   overlain
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
xPrim = np.arange(-rPrim, rPrim, 0.0001) 
yPrim = np.sqrt(np.square(rPrim) - np.square(xPrim)) 




#  convert Fracture pole plunges and plunge directions to cartesian coords    
dp = np.sqrt(1 - np.sin(fracture_poles_rad[:,0])) 
newTrend = np.remainder(fracture_poles_rad[:,1], 2*np.pi) 
xFractures = dp * np.sin(newTrend)  
yFractures = dp * np.cos(newTrend) 





#   convert principal stress orientations into cartesian coords for equal
#   area
[pstress, dCp] = principalstress(stress_tensor, trend_s1_rad, plunge_s1_rad, trend_s3_rad) 
trendS2rad = pstress[1,1] 
plungeS2rad = pstress[1,2] 
plungeS3rad = pstress[2,2] 
[xS1, yS1] = stcoordline(trend_s1_rad, plunge_s1_rad, 1)   
[xS2, yS2] = stcoordline(trendS2rad, plungeS2rad, 1)   
[xS3, yS3] = stcoordline(trend_s3_rad, plungeS3rad, 1)



    
        

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
#%   for each plane in the input file
#%   new loop to calculate specifc values for the supplied poles 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
sigmaNFracture = np.zeros((n_fractures,1))  
tauFracture = np.zeros((n_fractures,1))  
for i in range(0, n_fractures):
    
    #   convert fracture pole to strike and dip 
    [strike, dip] = pole(fracture_poles_rad[i,1], fracture_poles_rad[i,0], 0)  
    
    #   calculate normal and shear stress on the plane 
    [stressFracture, dcStressFracture, R] = shearonplane(stress_tensor, trend_s1_rad, plunge_s1_rad, trend_s3_rad, strike, dip) 
    
    #   save normal and shear stresses for later calculation 
    sigmaNFracture[i] = stressFracture[0,0]  ### Test this, i commented it out
    tauFracture[i] = stressFracture[2,0]  ### Test this, i commented it out
    
#   end for each plane 

#   calculate normalised slip tendency
TsFractureFile = (tauFracture / sigmaNFracture) / TsMax 
        
#   calculate dilatation tendency
TdFractureFile = ( sorted_sigma[0] - sigmaNFracture ) / ( sorted_sigma[0] - sorted_sigma[2] )     

#   calculate fracture susceptibility
SfFractureFile = sigmaNFracture - ( tauFracture / mu_static ) 

#   calculate opening angle
OAFile = tauFracture / (Pf - sigmaNFracture) 
muOAfractureFile = np.arctan(OAFile) * (180/np.pi) 

#   write out text file of data values for the specific fracture poles
all_data_for_export = np.hstack((fracture_poles, TsFractureFile, TdFractureFile, SfFractureFile,tauFracture, sigmaNFracture, muOAfractureFile))

np.savetxt("fracture_poles.csv", 
    all_data_for_export, 
    delimiter=',', 
    fmt='%2.2f'+','+'%3.2f'+','+'%1.3f'+','+'%1.3f'+','+'%5.2f'+','+'%5.2f'+','+'%5.2f'+','+'%5.3f', 
    header='plunge'+','+'trend'+','+'Ts'+','+'Td'+','+'Sf'+','+'tau'+','+'sigmaN'+','+'muOA', 
    comments='')






lwPrim = 1 
sizePoleMarker = 15 
sizeStressMarker = 10 
ncontours = 20 


def plot_stereonet_figs(dataset, n_fractures, ncontours, dataset_name='dataset_name'):

    '''
    Function to plot stereonet of different data.

    data is the input data you want to plot.
    1. slip tendency, Ts
    2. Dilatation Tendency, Td
    3. fracture susceptibility, Sf
    4. opening angle, OA

    dataset name is a string name for plot label.

    n_fractures is the number of fracture surfaces.

    ncontours is the number of contours in the plot.
    '''

    fig, ax = plt.subplots(1, figsize=(6,6))

    sn = plt.contourf(xeqarea, yeqarea, dataset, ncontours)#, 'EdgeColor', 'none')  
 
    ax.plot(xPrim, yPrim, '-k', linewidth=lwPrim)  
    ax.plot(xPrim, -yPrim, '-k', linewidth=lwPrim) 
    ax.plot(xFractures, yFractures, '.r', markersize=sizePoleMarker) 
    ax.plot(xS1, yS1, 's', markersize=sizeStressMarker, markeredgecolor='k', markerfacecolor='w')  
    ax.plot(xS2, yS2, 'd', markersize=sizeStressMarker, markeredgecolor='k', markerfacecolor='w')  
    ax.plot(xS3, yS3, '^', markersize=sizeStressMarker, markeredgecolor='k', markerfacecolor='w')
    ax.set_xticks([])
    ax.set_yticks([])
 
    plt.gca().set_aspect('equal', adjustable='box')

    plt.xlim([-1.05, 1.05])  
    plt.ylim([-1.05, 1.05])

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
 
    # creating colorbar
    fig.colorbar(sn, cax = cax, orientation = "vertical")

    fig.axes[0].set_title('{}'.format(dataset_name) + ' n={}'.format(n_fractures))

    #fig.colorbar(sn)
    #cb_ax = fig.axes[1] #Grab the second axis of the plot which is the colorbar
    # cb_ax.set_ylabel('{}'.format(dataset_name), 
    #              rotation=-90, 
    #              fontsize=14,
    #              ha='center',
    #              va='center',
    #              labelpad=10
    #             )
    
    fig.savefig("./figures/{}.pdf".format(dataset_name), dpi=300, transparent=True)
    plt.show()




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
#%  plot Mohr diagrams, with fractures & contoured stability measures  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

sizePoleMarker = 10  


def mohr_plot(sigma1, sigma2, sigma3, sigmaN, tau, dataset, ncontours, n_fractures, dataset_name='dataset_name'):

    '''
    Function to plot Mohr circles, uses plotmohr_OA.py

    sigma1: principal stresses
    sigma2: ''
    sigma3: ''
    sigmaN: numpy array with normal stresses
    tau: numpy array with shear stresses
    dataset: the data you want to plot
    ncountours: number of contours
    n_fractures: number of fractures from your dataset
    dataset_name: string of the name of your dataset for plotting

    '''




    ## Calculating the mohr circle
    theta_mohr = np.arange(0, 2*np.pi, 2*np.pi/360)
    sin2theta_mohr = np.sin(2 * theta_mohr)
    cos2theta_mohr = np.cos(2 * theta_mohr)

    tau13_mohr = ((sigma1 - sigma3)/2)*sin2theta_mohr
    sigma13_mohr = (sigma1 + sigma3)/2 + ((sigma1-sigma3)/2)*cos2theta_mohr

    tau12_mohr=((sigma1 - sigma2)/2)*sin2theta_mohr
    sigma12_mohr=(sigma1 + sigma2)/2 + ((sigma1 - sigma2)/2)*cos2theta_mohr

    tau23_mohr = ((sigma2 - sigma3)/2)*sin2theta_mohr
    sigma23_mohr = (sigma2 + sigma3)/2 + ((sigma2 - sigma3)/2)*cos2theta_mohr

    xm = sigma13_mohr.max()
    ym = tau13_mohr.max()
    

    ## Plot the mohr circle
    fig, ax = plt.subplots(1, figsize=(6,6)) 

    plt.plot([0, sigmaN_mohr], [cohesion, tau_mohr], '-r', linewidth=1)  
     
    ax.contourf(sigmaN, tau, dataset, ncontours)
    ax.plot(sigma13_mohr, tau13_mohr, linewidth=1, color='k')
    ax.plot(sigma12_mohr, tau12_mohr, linewidth=1, color='k')
    ax.plot(sigma23_mohr, tau23_mohr, linewidth=1, color='k')  


    for f in range(0,n_fractures): #1:n_fractures
        ax.plot(sigmaNFracture[f], tauFracture[f], '.r', markersize=sizePoleMarker) 
    
    plt.gca().set_aspect('equal', adjustable='box')

    plt.xlim([0, sigma1*1.05]) 
    plt.ylim([0, sigmad*0.75])
    ax.set_xlabel('Effective normal stress, MPa') 
    ax.set_ylabel('Shear stress, MPa') 
    plt.title('{}'.format(dataset_name)+' n='+'{}'.format(n_fractures)) 
    fig.savefig("./figures/{}.pdf".format(dataset_name), dpi=300, transparent=True) 

 
    plt.show()




plot_stereonet_figs(Ts, n_fractures, ncontours, dataset_name='Slip tendency, Ts') 
plot_stereonet_figs(Td, n_fractures, ncontours, dataset_name='Dilatation tendency, Td')
plot_stereonet_figs(Sf, n_fractures, ncontours, dataset_name='Fracture susceptibility, Sf')
plot_stereonet_figs(muOAfracture, n_fractures, ncontours, dataset_name='Opening angle, OA')



mohr_plot(sigma1, sigma2, sigma3, sigmaN, tau, Ts, ncontours, n_fractures, dataset_name='Slip tendency')
mohr_plot(sigma1, sigma2, sigma3, sigmaN, tau, Td, ncontours, n_fractures, dataset_name='Dilation tendency')
mohr_plot(sigma1, sigma2, sigma3, sigmaN, tau, Sf, ncontours, n_fractures, dataset_name='Fracture susceptibility')
mohr_plot(sigma1, sigma2, sigma3, sigmaN, tau, muOAfracture, ncontours, n_fractures, dataset_name='Opening angle')



#   2. Td 
# fig, ax = plt.subplot(figsize(6,6)) 
# set(gcf, 'PaperPositionMode', 'manual')  
# set(gcf, 'PaperUnits', 'inches') ; 
# set(gcf, 'PaperPosition', [ 0.25, 0.25, 5, 5])  

# plot([0, sigmaNMohr], [C0, tauMohr], '-r', 'LineWidth', 1)  
# #hold on ; 
# contourf(sigmaN, tau, Td, ncontours, 'EdgeColor', 'none')  
# plotMohr_OA(sigma1, sigma2, sigma3, 1, '-k')  
# cmocean('thermal') ; 
# for f in n_fractures: #= 1:n_fractures
#     plot(sigmaNFracture(f), tauFracture(f), '.r', 'MarkerSize', sizePoleMarker, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
# #hold off ; 
# xlim([0, sigma1*1.05])  
# ylim([0, sigmad*0.75]) 
# xlabel('Effective normal stress, MPa')  
# ylabel('Shear stress, MPa') 
# title('Dilatation tendency, T_d  n=', num2str(n_fractures)) 
# caxis([0, 1]) 
# cb = colorbar 
# cb.Location = 'SouthOutside' 
# cb.Label.String = 'Dilatation tendency'  
#print -r600 -dtiff 'FracTend_Td_mohr.tif' 
 
# #   3. Sf
# fig, ax = plt.subplot(figsize(6,6)) 
# set(gcf, 'PaperPositionMode', 'manual')  
# set(gcf, 'PaperUnits', 'inches')  
# set(gcf, 'PaperPosition', [ 0.25, 0.25, 5, 5]) 

# plot([0, sigmaNMohr], [C0, tauMohr], '-r', 'LineWidth', 1) 
# #hold on ; 
# contourf(sigmaN, tau, Sf, ncontours, 'EdgeColor', 'none')  
# plotMohr_OA(sigma1, sigma2, sigma3, 1, '-k')  
# colormap(flipud(cmocean('thermal'))) 
# for f in n_fractures: # = 1:n_fractures
#     plot(sigmaNFracture(f), tauFracture(f), '.r', 'MarkerSize', sizePoleMarker, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r') 

# #hold off ; 
# xlim([0, sigma1*1.05])  
# ylim([0, sigmad*0.75]) 
# xlabel('Effective normal stress, MPa') 
# ylabel('Shear stress, MPa') 
# title('Fracture susceptibility, S_f  n=', num2str(n_fractures)) 
# cb = colorbar  
# cb.Location = 'SouthOutside'  
# cb.Label.String = '\DeltaP_{f}, MPa'  
# #print -r600 -dtiff 'FracTend_Sf_mohr.tif' 

# #   opening angle Mohr diagram 
# fig, ax = plt.subplot(figsize(6,6)) 
# set(gcf, 'PaperPositionMode', 'manual') 
# set(gcf, 'PaperUnits', 'inches')
# set(gcf, 'PaperPosition', [ 0.25, 0.25, 5, 5])

# #   modified colourbar for opening angle scale 
# cbOA = cmocean('thermal', 9)  
# #cbOA(10, :) = [ 1, 1, 1 ]         #   white for OA < 0 
# cbOA = flipud(cbOA)               #   reverse it 

# plot([0, sigmaNMohr], [C0, tauMohr], '-r', 'LineWidth', 1)  
# #hold on ; 
# fill([0, Pf, Pf, 0], [0, 0, sigma1, sigma1], 'b', 'FaceAlpha', 0.2) 
# contourf(sigmaN, tau, muOAfracture, ncontours, 'EdgeColor', 'none') 
# plotMohr_OA(sigma1, sigma2, sigma3, 1, '-k')  
# for f in n_fractures: #= 1:n_fractures
#     plot(sigmaNFracture(f), tauFracture(f), '.r', 'MarkerSize', sizePoleMarker, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')  
# #hold off ; 
# colormap(cbOA) 
# xlim([0, sigma1*1.05])  
# ylim([0, sigmad*0.75]) 
# xlabel('Effective normal stress, MPa') 
# ylabel('Shear stress, MPa') 
# title('Opening angle \mu_a for P_f=', num2str(Pf), ' MPa  n=', num2str(n_fractures)) 
# caxis([-10, 90])
# cb = colorbar  
# cb.Location = 'SouthOutside'  
# cb.Label.String = 'Opening angle, \circ'  
# #print -r600 -dtiff 'FracTend_OA_mohr.tif'  
