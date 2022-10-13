
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable



def stereonet_plot(dataset, n_fractures, ncontours, dataset_name='dataset_name'):

    '''
    Function to plot stereonet of different data.

    data is the input data you want to plot.
    1. slip tendency, Ts
    2. Dilatation Tendency, Td
    3. fracture susceptibility, Sf
    4. opening angle, OA

    dataset name is a string name for plot label.

    n_fractures is the number of fracture surfaces.

    ncontours is the number of contours in the plot.
    '''

    fig, ax = plt.subplots(1, figsize=(6,6))

    sn = plt.contourf(xeqarea, yeqarea, dataset, ncontours)#, 'EdgeColor', 'none')  
 
    ax.plot(xPrim, yPrim, '-k', linewidth=lwPrim)  
    ax.plot(xPrim, -yPrim, '-k', linewidth=lwPrim) 
    ax.plot(xFractures, yFractures, '.r', markersize=sizePoleMarker) 
    ax.plot(xS1, yS1, 's', markersize=sizeStressMarker, markeredgecolor='k', markerfacecolor='w')  
    ax.plot(xS2, yS2, 'd', markersize=sizeStressMarker, markeredgecolor='k', markerfacecolor='w')  
    ax.plot(xS3, yS3, '^', markersize=sizeStressMarker, markeredgecolor='k', markerfacecolor='w')
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

    
    fig.savefig("./figures/{}.pdf".format(dataset_name), dpi=300, transparent=True)
    plt.show()




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
#%  plot Mohr diagrams, with fractures & contoured stability measures  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

sizePoleMarker = 10  


def mohr_plot(sigma1, sigma2, sigma3, sigmaN, tau, dataset, ncontours, n_fractures, dataset_name='dataset_name'):

    '''
    Function to plot Mohr circles, uses plotmohr_OA.py

    sigma1: principal stresses
    sigma2: ''
    sigma3: ''
    sigmaN: numpy array with normal stresses
    tau: numpy array with shear stresses
    dataset: the data you want to plot
    ncountours: number of contours
    n_fractures: number of fractures from your dataset
    dataset_name: string of the name of your dataset for plotting

    '''


    ## Calculating the mohr circle
    theta_mohr = np.arange(0, 2*np.pi, 2*np.pi/360)
    sin2theta_mohr = np.sin(2 * theta_mohr)
    cos2theta_mohr = np.cos(2 * theta_mohr)

    tau13_mohr = ((sigma1 - sigma3)/2)*sin2theta_mohr
    sigma13_mohr = (sigma1 + sigma3)/2 + ((sigma1-sigma3)/2)*cos2theta_mohr

    tau12_mohr=((sigma1 - sigma2)/2)*sin2theta_mohr
    sigma12_mohr=(sigma1 + sigma2)/2 + ((sigma1 - sigma2)/2)*cos2theta_mohr

    tau23_mohr = ((sigma2 - sigma3)/2)*sin2theta_mohr
    sigma23_mohr = (sigma2 + sigma3)/2 + ((sigma2 - sigma3)/2)*cos2theta_mohr

    xm = sigma13_mohr.max()
    ym = tau13_mohr.max()
    

    ## Plot the mohr circle
    fig, ax = plt.subplots(1, figsize=(6,6)) 

    plt.plot([0, sigmaN_mohr], [cohesion, tau_mohr], '-r', linewidth=1)  
     
    ax.contourf(sigmaN, tau, dataset, ncontours)
    ax.plot(sigma13_mohr, tau13_mohr, linewidth=1, color='k')
    ax.plot(sigma12_mohr, tau12_mohr, linewidth=1, color='k')
    ax.plot(sigma23_mohr, tau23_mohr, linewidth=1, color='k')  


    for f in range(0,n_fractures): 
        ax.plot(sigmaNFracture[f], tauFracture[f], '.r', markersize=sizePoleMarker) 
    
    plt.gca().set_aspect('equal', adjustable='box')

    plt.xlim([0, sigma1*1.05]) 
    plt.ylim([0, sigmad*0.75])
    ax.set_xlabel('Effective normal stress, MPa') 
    ax.set_ylabel('Shear stress, MPa') 
    plt.title('{}'.format(dataset_name)+' n='+'{}'.format(n_fractures)) 
    fig.savefig("./figures/{}.pdf".format(dataset_name), dpi=300, transparent=True) 

 
    plt.show()
