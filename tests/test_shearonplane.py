import numpy as np

from sga import shearonplane

def test_shearonplane():

	trend_s1 = np.radians(26) 
	plunge_s1 = np.radians(0)
	trend_s3 = np.radians(116)

	strike = np.radians(20)
	dip = np.radians(0)

	sigma1 = 50
	sigma2 = 20
	sigma3 = 5
	sorted_sigma = [sigma1, sigma2, sigma3]
	stress = np.diag(sorted_sigma)

	[TT,dCTT,R] = shearonplane(stress, trend_s1, plunge_s1, trend_s3, strike, dip)

test_shearonplane()

