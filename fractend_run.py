import os, sys, os.path
import numpy as np

from sga import stcoordline, pole, principalstress, shearonplane
from fractend import stress_state

import pdb


# fractures_file = r'C:\Users\bmelosh\Documents\Leapfrog\Raster_import\GM2\Polyline_MF2_vertices_out.csv'

# def ftp_run(fractures_file, output_name, sigma1, sigma2, sigma3, Pf, trend_s1, plunge_s1, mu_static, cohesion, sigmaN_mohr, ncontours, increment):

# 	## User inputs:

# 	#   read in poles to specific fractures; tab-delimited text file, 
# 	#    formatted as plunge then trend  
# 	#Utah_OA-Sills.txt'   # 'Utah_Thrusts_TS-RW. Utah_Sills_ArcRW. Utah_DefBands2RW. Utah_OA-Sills


# 	#fracture_poles = np.loadtxt(fractures_file, skiprows=1, delimiter=',')

# 	with open(fractures_file) as f:
# 	    #determining number of columns from the first line of text
# 	    n_cols = len(f.readline().split(","))
# 	fracture_poles = np.loadtxt(fractures_file, skiprows=1, delimiter="," , usecols=np.arange(n_cols-2, n_cols))

# 	np.savetxt(output_name+".csv", fracture_poles, delimiter=',', header='plunge'+','+'trend', comments='')



# if __name__=="__main__":
# 	fractures_file = sys.argv[1]
# 	ftp_run()



## User inputs:

#   read in poles to specific fractures; tab-delimited text file, 
#    formatted as plunge then trend  
output_filename = 'TEST2'
fractures_file = r'C:\Users\bmelosh\Documents\Leapfrog\Raster_import\GM2\Polyline_MF2_vertices_out.csv'#Utah_OA-Sills.txt'   # 'Utah_Thrusts_TS-RW. Utah_Sills_ArcRW. Utah_DefBands2RW. Utah_OA-Sills


#fracture_poles = np.loadtxt(fractures_file, skiprows=1, delimiter=',')

with open(fractures_file) as f:
    #determining number of columns from the first line of text
    n_cols = len(f.readline().split(","))
fracture_poles = np.loadtxt(fractures_file, skiprows=1, delimiter="," , usecols=np.arange(n_cols-2, n_cols))

np.savetxt(output_filename+".csv", fracture_poles, delimiter=',', header='plunge'+','+'trend', comments='')


#   read in stress magnitudes 
#   principal stresses in MPa 
sigma1 = 40       
sigma2 = 30     
sigma3 = 20 
# sigma1 = 100
# sigma2 = 48
# sigma3 = 29

Pf = 37 # pore fluid pressure in MPa, for fracture Opening Angle calculations


#   read in stress orientation 
#   e.g. for old case of SHmax azimuth of 135 
# normal fault system: Trend s1=0, Plunge s1=90 - use s3-trend to change orientation of stress field with s1 as the rotation axis
# thrust fault: Trend of s3 must be >90 from s1; Ps1=0
# strike-slip: Trend of s3 must be 90 from Trend of s1
# trend_s1 = 68 
# plunge_s1 =  3    
# trend_s3 = 265

trend_s1 = 26 #Stress orientations from Boyle and Zoback 2014
plunge_s1 = 0
trend_s3 = 116


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
#ss.mohr_plot(ss.Ts, 'Shear')



# Combine the data for import into leapfrog geo

fracture_poles_locations = np.loadtxt(fractures_file, skiprows=1, delimiter="," , usecols=np.arange(0, n_cols-2))

all_data = np.hstack((fracture_poles_locations, ss.all_data_for_export))

np.savetxt("FTP_results_alldata.csv", 
	all_data, 
	delimiter=',',
	fmt='%5.2f'+','+ '%6.2f'+',' + '%4.2f'+',' + '%2.2f'+',' + '%3.2f'+','+ '%1.0f'+',' +'%2.2f'+','+ '%3.2f'+','+'%3.2f'+','+'%1.3f'+','+'%1.3f'+','+'%5.2f'+','+'%5.2f'+','+'%5.2f'+','+'%5.3f', 
	header='x' + ',' + 'y' + ',' + 'z' + ',' + 'dip' + ',' + 'azimuth' + ',' + 'polarity' + ',' 'strike' + ',' + 'plunge'+','+'trend'+','+'Ts'+','+'Td'+','+'Sf'+','+'tau'+','+'sigmaN'+','+'muOA', 
	comments='')

