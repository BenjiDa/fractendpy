import numpy as np

from fractendpy.fractendpy_future import *
from plotting.figure_plotter import *


## User inputs:


#   read in poles to specific fractures; tab-delimited text file, 
#    formatted as plunge then trend 
fractures_file = 'Utah_OA-Sills.txt' ;  # 'Utah_Thrusts_TS-RW. Utah_Sills_ArcRW. Utah_DefBands2RW. Utah_OA-Sills

fracture_poles = np.loadtxt(fractures_file)
n_fractures = len(fracture_poles)
ncontours = 20

#   read in stress magnitudes 
#   principal stresses in MPa 
sigma1 = 43       
sigma2 = 38.86     
sigma3 = 25 

Pf = 37.24 # pore fluid pressure in MPa, for fracture Opening Angle calculations


#   read in stress orientation 
#   e.g. for old case of SHmax azimuth of 135 
# normal fault system: Trend s1=0, Plunge s1=90 - use s3-trend to change orientation of stress field with s1 as the rotation axis
# thrust fault: Trend of s3 must be >90 from s1; Ps1=0
# strike-slip: Trend of s3 must be 90 from Trend of s1
trend_s1 = 68 
plunge_s1 =  3    
trend_s3 = 265 

increment = 10 #how many degrees interval to calculate stresses, should be set to 10 for exploring data and 1 to calculate entire space. 

#   coefficient of friction & cohesion 
mu_static = 0.6 
cohesion = 0
sigmaN_mohr = 100
tau_mohr = cohesion + mu_static * sigmaN_mohr 

### End user inputs  



stress_tensor, sorted_sigma, sigmad = stress_tensor(sigma1, sigma2, sigma3)

sigmaN, tau = normal_and_shear_stress(stress_tensor, trend_s1, plunge_s1, trend_s3)

TsMax, Ts, Td, TD, Sf, OA, muOAfracture = tendencies(sigmaN, tau, sorted_sigma, mu_static, Pf)

#stress_ratio() ??

phiP, thetaP = azimuthal_variation(increment)

dp, xeqarea, yeqarea, rPrim, xPrim, yPrim = equal_area_projection(phiP, thetaP)

xFractures, yFractures = pole_to_cart(fracture_poles)

xS1, yS1, xS2, yS2, xS3, yS3 = stress_to_cart(stress_tensor, trend_s1, plunge_s1, trend_s3)

write_stress_to_file(fracture_poles, stress_tensor, sorted_sigma, TsMax, mu_static, Pf, trend_s1, plunge_s1, trend_s3)



###### Plotting ############
stereonet_plot(Ts, n_fractures, ncontours, dataset_name='Slip tendency, Ts') 
stereonet_plot(Td, n_fractures, ncontours, dataset_name='Dilatation tendency, Td')
stereonet_plot(Sf, n_fractures, ncontours, dataset_name='Fracture susceptibility, Sf')
stereonet_plot(muOAfracture, n_fractures, ncontours, dataset_name='Opening angle, OA')

mohr_plot(sigma1, sigma2, sigma3, sigmaN, tau, Ts, ncontours, n_fractures, dataset_name='Slip tendency')
mohr_plot(sigma1, sigma2, sigma3, sigmaN, tau, Td, ncontours, n_fractures, dataset_name='Dilation tendency')
mohr_plot(sigma1, sigma2, sigma3, sigmaN, tau, Sf, ncontours, n_fractures, dataset_name='Fracture susceptibility')
mohr_plot(sigma1, sigma2, sigma3, sigmaN, tau, muOAfracture, ncontours, n_fractures, dataset_name='Opening angle')
