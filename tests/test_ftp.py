import unittest
import numpy as np
from fractend_future_class import stress_state

import pdb

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



class test_stress_state(unittest.TestCase):

	## User inputs:

	fractures_file = 'Utah_OA-Sills.txt'   # 'Utah_Thrusts_TS-RW. Utah_Sills_ArcRW. Utah_DefBands2RW. Utah_OA-Sills
	fracture_poles = np.loadtxt(fractures_file)

	sigma1 = 43       
	sigma2 = 38.86     
	sigma3 = 25 

	Pf = 37.24 

	trend_s1 = 68 
	plunge_s1 =  3    
	trend_s3 = 265 

	mu_static = 0.6 
	cohesion = 0
	sigmaN_mohr = 100

	ncontours = 20 

	increment = 10 

	### End user inputs

	ss = stress_state(fracture_poles, sigma1, sigma2, sigma3, trend_s1, plunge_s1, trend_s3, Pf, mu_static, cohesion, sigmaN_mohr, increment, ncontours)
		


	def test_init(self):

		self.assertEqual(self.ss.fracture_poles.all(), self.fracture_poles.all())
		self.assertEqual(self.ss.sigma1, self.sigma1)
		self.assertEqual(self.ss.sigma2, self.sigma2)
		self.assertEqual(self.ss.sigma3, self.sigma3)
		self.assertEqual(self.ss.trend_s1, self.trend_s1)
		self.assertEqual(self.ss.plunge_s1, self.plunge_s1)
		self.assertEqual(self.ss.trend_s3, self.trend_s3)
		self.assertEqual(self.ss.Pf, self.Pf)
		self.assertEqual(self.ss.mu_static, self.mu_static)
		self.assertEqual(self.ss.cohesion, self.cohesion)
		self.assertEqual(self.ss.sigmaN_mohr, self.sigmaN_mohr)
		self.assertEqual(self.ss.increment, self.increment)
		self.assertEqual(self.ss.ncontours, self.ncontours)


	def test_stresstensor(self):

		sorted_sigma = [self.sigma1, self.sigma2, self.sigma3]
		self.assertEqual(self.ss.stresstensor().all(), np.diag(sorted_sigma).all())


	def test_normal_and_shear_stress(self):
		
		self.assertEqual(self.ss.trend_s1_rad, np.radians(self.trend_s1))
		self.assertEqual(self.ss.plunge_s1_rad, np.radians(self.plunge_s1))
		self.assertEqual(self.ss.trend_s3_rad, np.radians(self.trend_s3))

		self.assertEqual(self.ss.sigmaN.all(),  )

	def test_tendencies(self):
		pass

	def test_stress_ratios(self):
		pass

	def test_azimuthal_variation(self):
		pass

	def test_equal_area_projection(self):
		pass

	def test_pole_to_cart(self):
		pass

	def test_stress_to_cart(self):
		pass

	def test_write_stress_to_file(self):
		pass

	# def test_stereonet_plot(self, dataset, dataset_name='dataset_name'):
	# 	pass

	# def test_mohr_plot(self, dataset, dataset_name='dataset_name'):
	# 	pass



if __name__ == '__main__':
    unittest.main()