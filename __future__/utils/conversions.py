import numpy as np

'''
This file contains functions for fractendpy including zerotwopi, sphtocart, carttosph, and stcoordline
'''


def zerotwopi(a):

    '''
     zerotwopi constrains azimuth to lie between 0 and 2*pi radians
    
       b = ZeroTwoPi(a) returns azimuth b (from 0 to 2*pi)
       for input azimuth a (which may not be between 0 to 2*pi)
    
       NOTE: Azimuths a and b are input/output in radians 
    
    MATLAB script written by Nestor Cardozo for the book Structural 
    Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011. If you use
    this script, please cite this as "Cardozo in Allmendinger et al. (2011)"
    
    Converted to python by bmelosh Sept 27th 2022
    '''

    b=a

    twopi = 2.0*np.pi

    if b < 0.0:
        b = b + twopi
    elif b >= twopi:
        b = b - twopi

    return b

def sphtocart(trd,plg,k):
    

    '''
    sphtocart converts from spherical to cartesian coordinates 
    
       [cn,ce,cd] = SphToCart(trd,plg,k) returns the north (cn), 
       east (ce), and down (cd) direction cosines of a line.
    
       k is an integer to tell whether the trend and plunge of a line 
       (k = 0) or strike and dip of a plane in right hand rule 
       (k = 1) are being sent in the trd and plg slots. In this 
       last case, the direction cosines of the pole to the plane 
       are returned
    
       NOTE: Angles should be entered in radians 
    
    MATLAB script written by Nestor Cardozo for the book Structural 
    Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011. If you use
    this script, please cite this as "Cardozo in Allmendinger et al. (2011)"*/
    
    Converted to python by bmelosh Sept 28th 2022
    '''

    #If line (see Table 2.1)
    if k == 0:
        cd = np.sin(plg)
        ce = np.cos(plg) * np.sin(trd)
        cn = np.cos(plg) * np.cos(trd) 
    #Else pole to plane (see Table 2.1)
    elif k == 1:
        cd = np.cos(plg)
        ce = -1*np.sin(plg) * np.cos(trd)
        cn = np.sin(plg) * np.sin(trd)

    return [cn,ce,cd]


def carttosph(cn,ce,cd):

    '''
    carttosph converts from cartesian to spherical coordinates 
    
       [trd,plg] = carttosph(cn,ce,cd) returns the trend (trd)
       and plunge (plg) of a line for input north (cn), east (ce), 
       and down (cd) direction cosines
    
       NOTE: Trend and plunge are returned in radians
    
       CartToSph uses function ZeroTwoPi
    
    MATLAB script written by Nestor Cardozo for the book Structural 
    Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011. If you use
    this script, please cite this as "Cardozo in Allmendinger et al. (2011)"
    
    Converted to python by Ben Melosh Sept 28th 2022
    '''



    # Plunge (see Table 2.1)
    plg = np.arcsin(cd)

    # Trend
    # If north direction cosine is zero, trend is east or west
    # Choose which one by the sign of the east direction cosine
    if cn == 0.0: 
        if ce < 0.0:
            trd = 3.0/2.0*np.pi #% trend is west
        else:
            trd = np.pi/2.0 # trend is east

    # Else use Table 2.1
    else:
        trd = np.arctan(ce/cn) 
        if cn < 0.0:
            #Add pi 
            trd = trd+np.pi

        # Make sure trd is between 0 and 2*pi
        trd = zerotwopi(trd)

    return [trd,plg]



def stcoordline(trd, plg, sttype):

    '''
    stcoordline computes the coordinates of a line 
    in an equal angle or equal area stereonet of unit radius
     
        USE: [xp,yp] = stcoordline(trd,plg,sttype)
     
        trd = trend of line
        plg = plunge of line
        sttype = An integer indicating the type of stereonet. 0 for equal angle
                 and 1 for equal area
        xp and yp = Coordinates of the line in the stereonet plot
     
        NOTE: trend and plunge should be entered in radians
     
        stcoordlLine uses function zerotwopi
     
     MATLAB script written by Nestor Cardozo for the book Structural 
     Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011. If you use
     this script, please cite this as "Cardozo in Allmendinger et al. (2011)"
    
     Converted to python by Ben Melosh Sept 28th, 2022
     '''


    # Some constants
    piS4 = np.pi/4.0
    s2 = np.sqrt(2.0)
    plgS2 = plg/2.0

    # Take care of negative plunges
    if plg < 0.0:
        trd = zerotwopi(trd+np.pi)
        plg = -plg

    # Equal angle stereonet: From Equation 1.5 above
    # Also see Pollard and Fletcher (2005), eq.2.72
    if sttype == 0:
        xp = np.tan(piS4 - plgS2)*np.sin(trd)
        yp = np.tan(piS4 - plgS2)*np.cos(trd)
    # Equal area stereonet: From Equation 1.6 above
    # Also see Pollard and Fletcher (2005), eq.2.90
    elif sttype == 1:
        xp = s2*np.sin(piS4 - plgS2)*np.sin(trd)
        yp = s2*np.sin(piS4 - plgS2)*np.cos(trd)


    return [xp, yp]