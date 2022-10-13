import numpy as np

from utils import carttosph, sphtocart
from stress import principalstress, cauchy

def shearonplane(stress,tX1,pX1,tX3,strike,dip):

    '''
    shearonplane calculates the direction and magnitudes of the normal
    and shear tractions on an arbitrarily oriented plane
    
       USE: [TT,dCTT,R] = shearonplane(stress,tX1,pX1,tX3,strike,dip)
    
       stress = 3 x 3 stress tensor
       tX1 = trend of X1
       pX1 = plunge of X1
       tX3 = trend of X3
       strike = strike of plane
       dip = dip of plane
       TT = 3 x 3 matrix with the magnitude (column 1), trend (column 2) and 
           plunge (column 3) of: normal traction on the plane (row 1), 
           minimum shear traction (row 2), and maximum shear traction (row 3)
       dCTT = 3 x 3 matrix with the direction cosines of unit vectors parallel
             to: normal traction on the plane (row 1), minimum shear traction
             (row 2), and maximum shear traction (row 3)
       R = Stress ratio
    
       NOTE = Input stress tensor does not need to be along principal stress
              directions
              Plane orientation follows the right hand rule 
              Input/Output angles are in radians
    
       shearonplane uses functions principalstress, cauchy, carttosph, 
       and shptocart
    
    MATLAB script written by Nestor Cardozo for the book Structural 
    Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011. If you use
    this script, please cite this as "Cardozo in Allmendinger et al. (2011)"
    
    Converted to python by Ben Melosh Oct 3rd 2022
    '''

    # are there real or complex numbers
    w,v = np.linalg.eig(stress) 

    # Initialize TT and dCTT
    if np.iscomplex(v).any() == True:
        TT = np.zeros((3,3), dtype = np.complex_)
    else:
        TT = np.zeros((3,3))
    if np.iscomplex(v).any() == True:
        dCTT = np.zeros((3,3), dtype = np.complex_)
    else:
        dCTT = np.zeros((3,3))

    # Compute principal stresses and principal stress directions
    [pstress,dCp] = principalstress(stress,tX1,pX1,tX3)

    # Update stress vector so that it is along principal stress directions
    stress = np.zeros((3,3))
    stress = np.diag(pstress[:,0])

    #  New
    #  Calculate direction cosines of pole to plane
    if np.iscomplex(v).any() == True:
        p = np.zeros((1,3), dtype = np.complex_)
    else:
        p = np.zeros((1,3))
    [p[0][0],p[0][1],p[0][2]] = sphtocart(strike,dip,1)

    # Transform pole to plane to  principal stress coordinates
    if np.iscomplex(v).any() == True:
        pT = np.zeros((1,3), dtype = np.complex_)
    else:
        pT = np.zeros((1,3))
    # for i = 1:3
    #     for j = 1:3
    #         pT(i) = dCp(i,j)*p(j) + pT(i);
    pT = p @ dCp.transpose() + pT

    # Calculate the tractions in principal stress coordinates
    if np.iscomplex(v).any() == True:
        T = np.zeros((1,3), dtype = np.complex_)
    else:
        T = np.zeros((1,3))

    # Compute tractions using Cauchy's law
    #for i = 1:3
     #   for j = 1:3
     #       T(i) = stress(i,j)*pT(j) + T(i);

     # for i,j in range(1:3):
     #    T[i] = stress(i,j)*pT(j) + T(i)
    T = pT @ stress.transpose() + T


    # Find the B axis by the cross product of T cross pT and convert to
    # direction cosines (Eq 6.27)
    if np.iscomplex(v).any() == True:
        B = np.zeros((1,3), dtype = np.complex_)
    else:
        B = np.zeros((1,3))
    B[0][0] = T[0][1]*pT[0][2] - T[0][2]*pT[0][1]
    B[0][1] = T[0][2]*pT[0][0] - T[0][0]*pT[0][2]
    B[0][2] = T[0][0]*pT[0][1] - T[0][1]*pT[0][0]


    # Find the shear direction by the cross product of pT cross B. This will
    # give S in right handed coordinates (Eq. 6.27)
    if np.iscomplex(v).any() == True:
        S = np.zeros((1,3), dtype = np.complex_)
    else:
        S = np.zeros((1,3))
    S[0][0] = pT[0][1]*B[0][2] - pT[0][2]*B[0][1]
    S[0][1] = pT[0][2]*B[0][0] - pT[0][0]*B[0][2]
    S[0][2] = pT[0][0]*B[0][1] - pT[0][1]*B[0][0]

    # New: Convert B and S to unit vectors
    rB = np.sqrt(B[0][0]*B[0][0]+B[0][1]*B[0][1]+B[0][2]*B[0][2])
    rS = np.sqrt(S[0][0]*S[0][0]+S[0][1]*S[0][1]+S[0][2]*S[0][2])

    for i in range(0,3):
        B[0][i] = B[0][i]/rB
        S[0][i] = S[0][i]/rS

    # Now we can write the transformation matrix from principal stress
    # coordinates to plane coordinates (Eq. 6.28)
    if np.iscomplex(v).any() == True:
        aa = np.zeros((3,3), dtype = np.complex_)
    else:
        aa = np.zeros((3,3))
    aa[0] = [pT[0][0],pT[0][1],pT[0][2]]
    aa[1] = [B[0][0],B[0][1],B[0][2]]
    aa[2] = [S[0][0],S[0][1],S[0][2]]
    
    # Calculate stress ratio (Eq. 6.32)
    R = (stress[1][1] - stress[0][0])/(stress[2][2]-stress[0][0])

    # Calculate magnitude of normal and shear tractions (Eq. 6.31)
    for i in range(0,3):
        TT[i][0] = stress[0][0]*aa[0][0]*aa[i][0] + stress[1][1]*aa[0][1]*aa[i][1] + stress[2][2]*aa[0][2]*aa[i][2]

    # To get the orientation of the tractions in North-East-Down coordinates, we
    # need to do a vector transformation between principal stress and
    # North-East-Down coordinates. The transformation matrix are just the
    # direction cosines of the principal stresses in North-East-Down coordinates
    # (Eq. 6.29)

    # for i in range(0,3):
    #     for j in range(0,3):
    #         dCTT[0][i] = dCp[j][i]*pT[0][j] + dCTT[0][i]
    #         dCTT[1][i] = dCp[j][i]*B[0][j] + dCTT[1][i]
    #         dCTT[2][i] = dCp[j][i]*S[0][j] + dCTT[2][i]
    dCTT[0] = pT @ dCp + dCTT[0]
    dCTT[1] = B @ dCp + dCTT[1]
    dCTT[2] = S @ dCp + dCTT[2] 

    #Trend and plunge of traction on plane
    [TT[0][1],TT[0][2]] = carttosph(dCTT[0][0],dCTT[0][1],dCTT[0][2])
    #Trend and plunge of minimum shear direction
    [TT[1][1],TT[1][2]] = carttosph(dCTT[1][0],dCTT[1][1],dCTT[1][2])
    #Trend and plunge of maximum shear direction
    [TT[2][1],TT[2][2]] = carttosph(dCTT[2][0],dCTT[2][1],dCTT[2][2])


    return [TT,dCTT,R]