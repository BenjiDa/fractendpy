import numpy as np

from stress import principalstress

stress = np.array([[2,2,2],[2,2,2],[2,2,2]])
stress = np.array([[3, 1, 2],[2, 2, 2],[2, 4, 1]])

principalstress(stress, 1,1,1)