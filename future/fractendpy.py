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

#   read in poles to specific fractures; tab-delimited text file, 
#    formatted as plunge then trend 
print(' ') 
print('Reading input file of poles to planes...')   
fractures_file = 'Utah_OA-Sills.txt' ;  # 'Utah_Thrusts_TS-RW. Utah_Sills_ArcRW. Utah_DefBands2RW. Utah_OA-Sills

fracture_poles = np.loadtxt(fractures_file)
n_fractures = len(fracture_poles)
fracture_poles_rad = fracture_poles * np.pi / 180 
print('Read ', n_fractures, ' fracture poles')  




def stress_tensor(sigma1, sigma2, sigma3):#, Pf, trend_s1, plunge_s1, trend_s3):

    '''
    Returns 
    stress_tensor: 3x3 stress tensor with principal stresses
    sorted_sigma: 1X3 array with principal stresses
    sigmad: differential stress

    '''


    sorted_sigma = [sigma1, sigma2, sigma3] 
    sigmad = sigma1 - sigma3
    stress_tensor = np.diag(sorted_sigma) #   Generate stress tensor. note the implicit convention: x//s1, y//s2, z//s3
    print('Principal stresses ', sorted_sigma, ' in MPa')

    #trend_s1_rad = trend_s1 * np.pi / 180
    #plunge_s1_rad = plunge_s1 * np.pi / 180 
    #trend_s3_rad = trend_s3 * np.pi / 180


    return stress_tensor, sorted_sigma, sigmad




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





sorted_sigma = [sigma1, sigma2, sigma3]
sigmad = sigma1 - sigma3
stress_tensor = np.diag(sorted_sigma) #   Generate stress tensor. note the implicit convention: x//s1, y//s2, z//s3

#stress_tensor = xyxyx(43, 38.86, 25)       

 

  
print('Stress orientation:')                   
 




def normal_and_shear_stress(stress_tensor, trend_s1, plunge_s1, trend_s3):
    '''
    Calculate normal and shear stresses for all directions on all surfaces.

    Returns

    '''

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

    return sigmaN, tau
        

def tendencies(sigmaN, tau, sorted_sigma, mu_static):
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

    return Ts, Td, TD, Sf, OA, muOAfracture


def stress_ratios(sigma1, sigma2, sigma3):

    #   Calculate and display various stress ratios:
    Phi = ( sigma2 - sigma3 ) / ( sigma1 - sigma3 )  
    R = (sigma1 - sigma2) / (sigma1 - sigma3) 
    Rprime = ( Pf - sigma3 ) / ( sigma1 - sigma3 )  

    return Phi, R, Rprime


def azimuthal_variation(increment):

    #   plot azimuthal variation of tendencies, with poles to fractures
    #   overlain
    deltaP = increment * np.pi / 180 
    phiP = np.arange(np.pi/2, np.pi+deltaP, deltaP)
    phiP = phiP - np.pi/2  
    thetaP = np.arange(0, 2*np.pi+deltaP, deltaP)  
    [phiP, thetaP] = np.meshgrid(phiP, thetaP)

    return phiP, thetaP



def equal_area_projection(phiP, thetaP):

    #   equal area projection 
    dp = np.sqrt(1 - np.sin(phiP)) 
    xeqarea = dp * np.sin(thetaP)  
    yeqarea = dp * np.cos(thetaP)   
    rPrim = 1  
    xPrim = np.arange(-rPrim, rPrim, 0.0001) 
    yPrim = np.sqrt(np.square(rPrim) - np.square(xPrim)) 

    return dp, xeqarea, yeqarea, rPrim, xPrim, yPrim


def pole_to_cart(fracture_poles_rad):

    #  convert Fracture pole plunges and plunge directions to cartesian coords    
    dp = np.sqrt(1 - np.sin(fracture_poles_rad[:,0])) 
    newTrend = np.remainder(fracture_poles_rad[:,1], 2*np.pi) 
    xFractures = dp * np.sin(newTrend)  
    yFractures = dp * np.cos(newTrend) 

    return xFractures, yFractures



def stress_to_cart(stress_tensor, trend_s1, plunge_s1, trend_s3):

    '''
    '''


    # convert stress orientations into radians
    trend_s1_rad = trend_s1 * np.pi / 180
    plunge_s1_rad = plunge_s1 * np.pi / 180 
    trend_s3_rad = trend_s3 * np.pi / 180

    #   convert principal stress orientations into cartesian coords for equal
    #   area
    [pstress, dCp] = principalstress(stress_tensor, trend_s1_rad, plunge_s1_rad, trend_s3_rad) 
    trendS2rad = pstress[1,1] 
    plungeS2rad = pstress[1,2] 
    plungeS3rad = pstress[2,2] 
    [xS1, yS1] = stcoordline(trend_s1_rad, plunge_s1_rad, 1)   
    [xS2, yS2] = stcoordline(trendS2rad, plungeS2rad, 1)   
    [xS3, yS3] = stcoordline(trend_s3_rad, plungeS3rad, 1)

    return    xS1, yS1, xS2, yS2, xS3, yS3



def write_stress_to_file():
    
    '''
    Calculates stress and shear values for each supplied pole, returns .csv file
    '''

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
    SfFractureFile = sigmaNFracture - ( tauFracture / muStatic ) 

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