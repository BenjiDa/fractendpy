import numpy as np

from sga import stcoordline, pole, principalstress, shearonplane
from fractend import stress_state

## User inputs:

#   read in poles to specific fractures; tab-delimited text file, 
#    formatted as plunge then trend  
fractures_file = 'Utah_OA-Sills.txt'   # 'Utah_Thrusts_TS-RW. Utah_Sills_ArcRW. Utah_DefBands2RW. Utah_OA-Sills

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

ncontours = 20 #number of contours

increment = 10 #increment to do calculations by, with a range of 0 to 360 and 90 to 180 for 3d space. 
# value should be set to 10 for testing and 1 for final run.

### End user inputs

 

## Run program

ss = stress_state(fracture_poles, sigma1, sigma2, sigma3, trend_s1, plunge_s1, trend_s3, Pf, mu_static, cohesion, sigmaN_mohr, increment, ncontours)

ss.stereonet_plot(ss.Ts, 'shear')
ss.mohr_plot(ss.Ts, 'Shear')

