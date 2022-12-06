# fractendpy
Python version of [FracTend.m](https://github.com/DaveHealy-github/FracTend)

Plotting functions are working but tests are not complete, work in progress. Outputs include figures and .csv file of data.

## To run the code: 
Clone the repository, navigate to the directory, and run "python fractend_run.py". This will run the code for the data file that is included 'Utah_OA-Sills'.
To run on your own data open fractend_run.py in an editor and enter your data file path, stress magnitudes, orientations,
pore fluid pressure, cohesion, etc. Then save and run "python fractend_run.py" in the command line while in the correct directory.


## Examples:

```
fractures_file = r'C:\Users\path\to\your\data.csv'   # Must be two columns of plunge and trend of pole to fault plane
fracture_poles = np.loadtxt(fractures_file, skiprows=1, delimiter=',')

# Example user imputs
#   read in stress magnitudes 
#   principal stresses in MPa 
sigma1 = 50       
sigma2 = 20     
sigma3 = 5 

Pf = 37 # pore fluid pressure in MPa, for fracture Opening Angle calculations

#   read in stress orientation 
#   e.g. for old case of SHmax azimuth of 135 
# normal fault system: Trend s1=0, Plunge s1=90 - use s3-trend to change orientation of stress field with s1 as the rotation axis
# thrust fault: Trend of s3 must be >90 from s1; Ps1=0
# strike-slip: Trend of s3 must be 90 from Trend of s1
trend_s1 = 26
plunge_s1 = 0
trend_s3 = 116

#   coefficient of friction & cohesion 
mu_static = 0.6 
cohesion = 10
sigmaN_mohr = 100

ncontours = 20 #number of contours, used for plotting

increment = 10 #increment to do calculations by, with a range of 0 to 360 and 90 to 180 for 3d space. 
# value should be set to 10 for testing and 1 for final run.
# End user inputs

# Run program

ss = stress_state(fracture_poles, sigma1, sigma2, sigma3, trend_s1, plunge_s1, trend_s3, Pf, mu_static, cohesion, sigmaN_mohr, increment, ncontours)

ss.stereonet_plot(ss.tau, 'Shear stress (MPa)')
ss.stereonet_plot(ss.sigmaN, 'Normal stress (MPa)')
ss.mohr_plot(ss.sigmaN, 'Normal stress') # plotting functions for mohr circles don't work well yet.
```

Results include a data table with all results for each fault surface and plots.


## To do: 

Finish tests for structural geology algorithms (sga.py) code from Cardozo and Allmendinger et al. (2012)
Finish tests for fractend



