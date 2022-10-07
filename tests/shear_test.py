import numpy as np

from shear import shearonplane

stress = np.array([[2,2,2],[2,2,2],[2,2,2]]) #numbers very small, does not match exactly with matlab results
stress = np.array([[3, 1, 2],[2, 2, 2],[2, 4, 1]]) #matches matlab results, seems to work fine
stress = np.array([[1, 1, 4],[4, 3, 6],[1, 5, 2]]) #complex
[TT,dCTT,R] = shearonplane(stress, 1,1,1,1,1)

print('TT: ',TT)

print('dCTT: ',dCTT)

print('R: ',R)