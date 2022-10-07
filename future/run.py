from fractendpy import fractendpy




#   read in poles to specific fractures; tab-delimited text file, 
#    formatted as plunge then trend 
fractures_file = 'Utah_OA-Sills.txt' 
fracture_poles = np.loadtxt(fractures_file)
n_fractures = len(fracture_poles)
fracture_poles_rad = fracture_poles * np.pi / 180 


## User inputs:

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

## End user inputs




stress_tensor, sorted_sigma, sigmad = stress_tensor(sigma1, sigma2, sigma3)

sigmaN, tau = normal_and_shear_stress(stress_tensor, trend_s1, plunge_s1, trend_s3)

Ts, Td, TD, Sf, OA, muOAfracture = tendencies(sigmaN, tau, sorted_sigma, mu_static)