import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from sga import stcoordline, pole, principalstress, shearonplane

'''
    fractendpy is the python version of FracTend by David Healy & Tara Stephens (July 2018), it runs on a base of code from
    Allmendinger et al., 2012 "Structural geology algorithms".  Additional equations and code from Morris et al., 1996 Geology; 
    Streit & Hillis, 2004 Energy; Jolly & Sanderson, 1997 Journal of Structural Geology. 

    Requires numpy and matplotlib

    Oct 2022
    Ben Melosh
'''




class stress_state():

    def __init__(self, fracture_poles, sigma1, sigma2, sigma3, trend_s1, plunge_s1, trend_s3, Pf, mu_static, cohesion, sigmaN_mohr, increment, ncontours):
        # User inputs
        self.fracture_poles = fracture_poles
        self.sigma1 = sigma1
        self.sigma2 = sigma2
        self.sigma3 = sigma3
        self.trend_s1 = trend_s1
        self.plunge_s1 = plunge_s1
        self.trend_s3 = trend_s3
        self.Pf = Pf
        self.mu_static = mu_static
        self.cohesion = cohesion
        self.sigmaN_mohr = sigmaN_mohr
        self.increment = increment
        self.ncontours = ncontours

        # Run
        self.stresstensor()
        self.normal_and_shear_stress()
        self.tendencies()
        self.stress_ratios()
        self.azimuthal_variation()
        self.equal_area_projection()
        self.pole_to_cart()
        self.stress_to_cart()
        self.write_stress_to_file()


    def stresstensor(self):

        '''
        Returns 
        stress_tensor: 3x3 stress tensor with principal stresses
        sorted_sigma: 1X3 array with principal stresses
        sigmad: differential stress

        '''


        self.sorted_sigma = [self.sigma1, self.sigma2, self.sigma3]
        self.sigmad = self.sigma1 - self.sigma3
        self.stress_tensor = np.diag(self.sorted_sigma)

        return self.stress_tensor
        


    def normal_and_shear_stress(self):
        '''
        Calculate normal and shear stresses for all directions on all surfaces.

        Returns

        '''

        # convert stress orientations into radians
        self.trend_s1_rad = np.radians(self.trend_s1)# * np.pi / 180
        self.plunge_s1_rad = np.radians(self.plunge_s1)# * np.pi / 180 
        self.trend_s3_rad = np.radians(self.trend_s3)# * np.pi / 180


        #   define a 3d space
        # for all directions in 3-space, calculate normal and shear stresses 
        #   on all 3d surfaces    
        phi_index = 0  
        theta_index = 0  
        #increment = 10
        phi_min = 90  
        phi_max = 180  
        theta_min = 0  
        theta_max = 360  
        phi_n = (phi_max - phi_min)//self.increment + 1  #divide to make int
        theta_n = (theta_max - theta_min)//self.increment + 1 
        self.sigmaN = np.zeros((theta_n, phi_n))  
        self.tau = np.zeros((theta_n, phi_n)) 


        for idp, phi in enumerate(range(90,181,self.increment)):

            phi_rad = (phi - 90) * np.pi/180

            for idt, theta in enumerate(range(0,361,self.increment)):

                theta_rad = theta * np.pi / 180

                #   convert pole to strike and dip 
                [strike, dip] = pole(theta_rad, phi_rad, 0) 

                #   calculate normal and shear stress on the plane
                [stress_fracture, dc_stress_fracture, R] = shearonplane(self.stress_tensor, self.trend_s1_rad, self.plunge_s1_rad, self.trend_s3_rad, strike, dip) 

                np.nan_to_num(stress_fracture, copy=False, nan=0.0) # This line converts nan values into zero, nan values are generated when
                #the dip of the surface is parallel to sigma 1 and sigma 3

                #   save normal and shear stresses for later calculation
                self.sigmaN[idt][idp] = stress_fracture[0][0]#(1,1) 
                self.tau[idt][idp] = stress_fracture[2][0]#(3,1)

        # import pdb
        # pdb.set_trace()
        #return sigmaN, tau
        
        

    def tendencies(self):#sigmaN, tau, sorted_sigma, mu_static, Pf):
        #   calculate tendencies - slip, dilatation and frac. suscep
        
        #   calculate normalised slip tendency (Morris et al., 1996)
        self.TsMax = (self.tau / self.sigmaN ).max() 
        self.Ts = ( self.tau / self.sigmaN ) / self.TsMax  

        #   calculate dilatation tendency (Ferril et al., 1999)
        self.Td =  ( self.sorted_sigma[0] - self.sigmaN ) / ( self.sorted_sigma[0] - self.sorted_sigma[2] )   

        #   calculate shear stress/ dilation displacement ratio 
        #   from Delaney et al.(1988)
        self.TD = self.tau / (self.sorted_sigma[2] + self.Pf) 

        #   calculate fracture susceptibility
        self.Sf = self.sigmaN - ( self.tau / self.mu_static )     

        #   calculate muOA (opening angle), Jolly & Sanderson, 1997
        OA = self.tau / (self.Pf - self.sigmaN)
        self.muOAfracture = np.arctan(OA) * (180/np.pi)

        #return TsMax, Ts, Td, TD, Sf, OA, muOAfracture


    def stress_ratios(self):

        #   Calculate and display various stress ratios:
        Phi = ( self.sigma2 - self.sigma3 ) / ( self.sigma1 - self.sigma3 )  
        R = (self.sigma1 - self.sigma2) / (self.sigma1 - self.sigma3) 
        Rprime = ( self.Pf - self.sigma3 ) / ( self.sigma1 - self.sigma3 )  
        print('Printing some stress ratios')
        print('---------------------------')
        print('Phi: {:1.2f}'.format(Phi))
        print('R: {:1.2f}'.format(R))
        print('Rprime: {:1.2f}'.format(Rprime))
        print('---------------------------')

        #return Phi, R, Rprime


    def azimuthal_variation(self):

        #   plot azimuthal variation of tendencies, with poles to fractures
        #   overlain
        deltaP = self.increment * np.pi / 180 
        phiP = np.arange(np.pi/2, np.pi+deltaP, deltaP)
        phiP = phiP - np.pi/2  
        thetaP = np.arange(0, 2*np.pi+deltaP, deltaP)  
        [self.phiP, self.thetaP] = np.meshgrid(phiP, thetaP)

        #return phiP, thetaP



    def equal_area_projection(self):

        #   equal area projection 
        dp = np.sqrt(1 - np.sin(self.phiP)) 
        self.xeqarea = dp * np.sin(self.thetaP)  
        self.yeqarea = dp * np.cos(self.thetaP)   
        self.rPrim = 1  
        self.xPrim = np.arange(self.rPrim*-1, self.rPrim, 0.0001) 
        self.yPrim = np.sqrt(np.square(self.rPrim) - np.square(self.xPrim)) 

        #return dp, xeqarea, yeqarea, rPrim, xPrim, yPrim


    def pole_to_cart(self):


        fracture_poles_rad = np.radians(self.fracture_poles)# * np.pi / 180

        #  convert Fracture pole plunges and plunge directions to cartesian coords    
        dp = np.sqrt(1 - np.sin(fracture_poles_rad[:,0])) 
        newTrend = np.remainder(fracture_poles_rad[:,1], 2*np.pi) 
        self.xFractures = dp * np.sin(newTrend)  
        self.yFractures = dp * np.cos(newTrend) 

        #return xFractures, yFractures


    def stress_to_cart(self):

        '''
        '''


        # convert stress orientations into radians
        trend_s1_rad = np.radians(self.trend_s1) 
        plunge_s1_rad = np.radians(self.plunge_s1)  
        trend_s3_rad = np.radians(self.trend_s3) 

        #   convert principal stress orientations into cartesian coords for equal
        #   area
        [pstress, dCp] = principalstress(self.stress_tensor, trend_s1_rad, plunge_s1_rad, trend_s3_rad) 
        trendS2rad = pstress[1,1] 
        plungeS2rad = pstress[1,2] 
        plungeS3rad = pstress[2,2] 
        [self.xS1, self.yS1] = stcoordline(trend_s1_rad, plunge_s1_rad, 1)   
        [self.xS2, self.yS2] = stcoordline(trendS2rad, plungeS2rad, 1)   
        [self.xS3, self.yS3] = stcoordline(trend_s3_rad, plungeS3rad, 1)

        # import pdb
        # pdb.set_trace()

        #return    xS1, yS1, xS2, yS2, xS3, yS3



    def write_stress_to_file(self):
        
        '''
        Calculates stress and shear values for each supplied pole, returns .csv file
        '''

        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        #%   for each plane in the input file
        #%   new loop to calculate specifc values for the supplied poles 
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

        fracture_poles_rad = np.radians(self.fracture_poles) #* np.pi / 180
        n_fractures = len(self.fracture_poles)

        # convert stress orientations into radians
        trend_s1_rad = np.radians(self.trend_s1) # * np.pi / 180
        plunge_s1_rad = np.radians(self.plunge_s1) #* np.pi / 180 
        trend_s3_rad = np.radians(self.trend_s3) #* np.pi / 180


        self.sigmaNFracture = np.zeros((n_fractures,1))  
        self.tauFracture = np.zeros((n_fractures,1))  
        for i in range(0, n_fractures):
            
            #   convert fracture pole to strike and dip 
            [strike, dip] = pole(fracture_poles_rad[i,1], fracture_poles_rad[i,0], 0)  
            
            #   calculate normal and shear stress on the plane 
            [stressFracture, dcStressFracture, R] = shearonplane(self.stress_tensor, trend_s1_rad, plunge_s1_rad, trend_s3_rad, strike, dip) 
            
            #   save normal and shear stresses for later calculation 
            self.sigmaNFracture[i] = stressFracture[0,0]  
            self.tauFracture[i] = stressFracture[2,0]  
        #   end for each plane 

        #   calculate normalised slip tendency
        TsFractureFile = (self.tauFracture / self.sigmaNFracture) / self.TsMax 
                
        #   calculate dilatation tendency
        TdFractureFile = ( self.sorted_sigma[0] - self.sigmaNFracture ) / ( self.sorted_sigma[0] - self.sorted_sigma[2] )     

        #   calculate fracture susceptibility
        SfFractureFile = self.sigmaNFracture - ( self.tauFracture / self.mu_static ) 

        #   calculate opening angle
        OAFile = self.tauFracture / (self.Pf - self.sigmaNFracture) 
        muOAfractureFile = np.arctan(OAFile) * (180/np.pi) 

        #   write out text file of data values for the specific fracture poles
        self.all_data_for_export = np.hstack((self.fracture_poles, TsFractureFile, TdFractureFile, SfFractureFile, self.tauFracture, self.sigmaNFracture, muOAfractureFile))

        np.savetxt("ftp_results.csv", 
            self.all_data_for_export, 
            delimiter=',', 
            fmt='%2.2f'+','+'%3.2f'+','+'%1.3f'+','+'%1.3f'+','+'%5.2f'+','+'%5.2f'+','+'%5.2f'+','+'%5.3f', 
            header='plunge'+','+'trend'+','+'Ts'+','+'Td'+','+'Sf'+','+'tau'+','+'sigmaN'+','+'muOA', 
            comments='')

        return self.all_data_for_export


    def stereonet_plot(self, dataset, dataset_name='dataset_name'):#stress_tensor, stress_orientations, increment, fracture_poles, ncontours, '):

        '''
        Function to plot stereonet of different tendency data.

        dataset is the input data you want to plot.
            slip tendency, Ts
            Dilatation Tendency, Td
            fracture susceptibility, Sf
            opening angle, OA
        stress_tensor: 3x3 array of principle stress magnitudes
        stress_orientations: list of trends and plunge of stress orientations
        increment: value defining the resolution to calculate and plot at.
        fracture_poles: plunge, trend values of poles to fracture planes from initial input.
        ncontours: the number of contours in the plot.
        dataset_name: a string name for plot label.

        Ben Melosh Oct 2022
        '''

        n_fractures = len(self.fracture_poles) # Calculate number of fractures
        # calculate plotting area and tendency variations

        fig, ax = plt.subplots(1, figsize=(6,6))

        sn = plt.contourf(self.xeqarea, self.yeqarea, dataset, self.ncontours)#, 'EdgeColor', 'none')  
        
        ax.plot(self.xPrim, self.yPrim, '-k', linewidth=1)  
        ax.plot(self.xPrim, self.yPrim*-1, '-k', linewidth=1) 
        ax.plot(self.xFractures, self.yFractures, '.r',markeredgecolor='k', markersize=15) 
        ax.plot(self.xS1, self.yS1, 's', markersize=10, markeredgecolor='k', markerfacecolor='w')
        ax.annotate(text='\u03C3'+'1', xy=(self.xS1+0.05, self.yS1+0.05))  
        ax.plot(self.xS2, self.yS2, 'd', markersize=10, markeredgecolor='k', markerfacecolor='w') 
        ax.annotate(text='\u03C3'+'2', xy=(self.xS2+0.05, self.yS2+0.05))  
        ax.plot(self.xS3, self.yS3, '^', markersize=10, markeredgecolor='k', markerfacecolor='w')
        ax.annotate(text='\u03C3'+'3', xy=(self.xS3+0.05, self.yS3+0.05)) 
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
        
        fig.savefig("./{}_stereo.pdf".format(dataset_name), dpi=300, transparent=True)
        plt.show()


    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    #%  plot Mohr diagrams, with fractures & contoured stability measures  
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


    def mohr_plot(self, dataset, dataset_name='dataset_name'):  

        '''
        Function to plot Mohr circles, uses plotmohr_OA.py
        
        dataset: the tendency data you want to plot, 
            slip tendency, Ts
            Dilatation Tendency, Td
            fracture susceptibility, Sf
            opening angle, OA
        stress_tensor: 3x3 tensor of the principal stresses
        sigmaN: numpy array with normal stresses
        tau: numpy array with shear stresses
        ncountours: number of contours
        fracture_poles: number of fractures from your dataset
        dataset_name: string of the name of your dataset for plotting

        '''


        n_fractures = len(self.fracture_poles)

        #sigmaN_mohr = 100
        tau_mohr = self.cohesion + self.mu_static * self.sigmaN_mohr 

        ## Calculating the mohr circle
        theta_mohr = np.arange(0, 2*np.pi, 2*np.pi/360)
        sin2theta_mohr = np.sin(2 * theta_mohr)
        cos2theta_mohr = np.cos(2 * theta_mohr)

        tau13_mohr = ((self.sigma1 - self.sigma3)/2)*sin2theta_mohr
        sigma13_mohr = (self.sigma1 + self.sigma3)/2 + ((self.sigma1-self.sigma3)/2)*cos2theta_mohr

        tau12_mohr=((self.sigma1 - self.sigma2)/2)*sin2theta_mohr
        sigma12_mohr=(self.sigma1 + self.sigma2)/2 + ((self.sigma1 - self.sigma2)/2)*cos2theta_mohr

        tau23_mohr = ((self.sigma2 - self.sigma3)/2)*sin2theta_mohr
        sigma23_mohr = (self.sigma2 + self.sigma3)/2 + ((self.sigma2 - self.sigma3)/2)*cos2theta_mohr

        xm = sigma13_mohr.max()
        ym = tau13_mohr.max()
        

        ## Plot the mohr circle
        fig, ax = plt.subplots(1, figsize=(6,6)) 

        plt.plot([0, self.sigmaN_mohr], [self.cohesion, tau_mohr], '-r', linewidth=1)  
         
        ax.contourf(self.sigmaN, self.tau, dataset, self.ncontours)
        ax.plot(sigma13_mohr, tau13_mohr, linewidth=1, color='k')
        ax.plot(sigma12_mohr, tau12_mohr, linewidth=1, color='k')
        ax.plot(sigma23_mohr, tau23_mohr, linewidth=1, color='k')  

        for f in range(0,n_fractures): #1:n_fractures
            ax.plot(self.sigmaNFracture[f], self.tauFracture[f], '.r', markersize=10) 
        
        plt.gca().set_aspect('equal', adjustable='box')

        plt.xlim([0, self.sigma1*1.05]) 
        plt.ylim([0, self.sigmad*0.75])
        ax.set_xlabel('Effective normal stress, MPa') 
        ax.set_ylabel('Shear stress, MPa') 
        plt.title('{}'.format(dataset_name)+' n='+'{}'.format(n_fractures)) 
        fig.savefig("./{}_mohr.pdf".format(dataset_name), dpi=300, transparent=True) 

     
        plt.show()
