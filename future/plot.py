

lwPrim = 1 
sizePoleMarker = 15 
sizeStressMarker = 10 
ncontours = 20 


def plot_stereonet_figs(dataset, n_fractures, ncontours, dataset_name='dataset_name'):

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

    xS1, yS1, xS2, yS2, xS3, yS3 = stress_to_cart(stress_tensor, trend_s1, plunge_s1, trend_s3)
  

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

    #fig.colorbar(sn)
    #cb_ax = fig.axes[1] #Grab the second axis of the plot which is the colorbar
    # cb_ax.set_ylabel('{}'.format(dataset_name), 
    #              rotation=-90, 
    #              fontsize=14,
    #              ha='center',
    #              va='center',
    #              labelpad=10
    #             )
    
    fig.savefig("{}.pdf".format(dataset_name), dpi=300, transparent=True)
    plt.show()






#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
#%  plot Mohr diagrams, with fractures & contoured stability measures  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

sizePoleMarker = 10  

#   coefficient of friction & cohesion 
mu_static = 0.6 
cohesion = 0 
sigmaN_mohr = 100
tau_mohr = cohesion + mu_static * sigmaN_mohr 


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
    sigma12_mohr=(sigma2 + sigma3)/2 + ((sigma2 - sigma3)/2)*cos2theta_mohr

    tau23_mohr = ((sigma2 - sigma3)/2)*sin2theta_mohr
    sigma23_mohr = (sigma2 + sigma3)/2 + ((sigma2 - sigma3)/2)*cos2theta_mohr

    xm = sigma13_mohr.max()
    ym = tau13_mohr.max()
    

    ## Plot the mohr circle
    fig, ax = plt.subplots(1, figsize=(6,6)) 

    plt.plot([0, sigmaN_mohr], [cohesion, tau_mohr], '-r', linewidth=1)  
     
    ax.contourf(sigmaN, tau, dataset, ncontours)
    ax.plot(sigma13_mohr, tau13_mohr, linewidth=1)
    ax.plot(sigma12_mohr, tau12_mohr, linewidth=1)
    ax.plot(sigma23_mohr, tau23_mohr, linewidth=1)  


    for f in range(0,n_fractures): #1:n_fractures
        ax.plot(sigmaNFracture[f], tauFracture[f], '.r', markersize=sizePoleMarker) 
    
    plt.gca().set_aspect('equal', adjustable='box')

    plt.xlim([0, sigma1*1.05]) 
    plt.ylim([0, sigmad*0.75])
    ax.set_xlabel('Effective normal stress, MPa') 
    ax.set_ylabel('Shear stress, MPa') 
    plt.title('{}'.format(dataset_name)+' n='+'{}'.format(n_fractures))  

 
    plt.show()


mohr_plot(sigma1, sigma2, sigma3, sigmaN, tau, Ts, 1000, n_fractures)



#   2. Td 
fig, ax = plt.subplot(figsize(6,6)) 
set(gcf, 'PaperPositionMode', 'manual')  
set(gcf, 'PaperUnits', 'inches') ; 
set(gcf, 'PaperPosition', [ 0.25, 0.25, 5, 5])  

plot([0, sigmaNMohr], [C0, tauMohr], '-r', 'LineWidth', 1)  
#hold on ; 
contourf(sigmaN, tau, Td, ncontours, 'EdgeColor', 'none')  
plotMohr_OA(sigma1, sigma2, sigma3, 1, '-k')  
cmocean('thermal') ; 
for f in n_fractures: #= 1:n_fractures
    plot(sigmaNFracture(f), tauFracture(f), '.r', 'MarkerSize', sizePoleMarker, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
#hold off ; 
xlim([0, sigma1*1.05])  
ylim([0, sigmad*0.75]) 
xlabel('Effective normal stress, MPa')  
ylabel('Shear stress, MPa') 
title('Dilatation tendency, T_d  n=', num2str(n_fractures)) 
caxis([0, 1]) 
cb = colorbar 
cb.Location = 'SouthOutside' 
cb.Label.String = 'Dilatation tendency'  
#print -r600 -dtiff 'FracTend_Td_mohr.tif' 
 
#   3. Sf
fig, ax = plt.subplot(figsize(6,6)) 
set(gcf, 'PaperPositionMode', 'manual')  
set(gcf, 'PaperUnits', 'inches')  
set(gcf, 'PaperPosition', [ 0.25, 0.25, 5, 5]) 

plot([0, sigmaNMohr], [C0, tauMohr], '-r', 'LineWidth', 1) 
#hold on ; 
contourf(sigmaN, tau, Sf, ncontours, 'EdgeColor', 'none')  
plotMohr_OA(sigma1, sigma2, sigma3, 1, '-k')  
colormap(flipud(cmocean('thermal'))) 
for f in n_fractures: # = 1:n_fractures
    plot(sigmaNFracture(f), tauFracture(f), '.r', 'MarkerSize', sizePoleMarker, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r') 

#hold off ; 
xlim([0, sigma1*1.05])  
ylim([0, sigmad*0.75]) 
xlabel('Effective normal stress, MPa') 
ylabel('Shear stress, MPa') 
title('Fracture susceptibility, S_f  n=', num2str(n_fractures)) 
cb = colorbar  
cb.Location = 'SouthOutside'  
cb.Label.String = '\DeltaP_{f}, MPa'  
#print -r600 -dtiff 'FracTend_Sf_mohr.tif' 

#   opening angle Mohr diagram 
fig, ax = plt.subplot(figsize(6,6)) 
set(gcf, 'PaperPositionMode', 'manual') 
set(gcf, 'PaperUnits', 'inches')
set(gcf, 'PaperPosition', [ 0.25, 0.25, 5, 5])

#   modified colourbar for opening angle scale 
cbOA = cmocean('thermal', 9)  
#cbOA(10, :) = [ 1, 1, 1 ]         #   white for OA < 0 
cbOA = flipud(cbOA)               #   reverse it 

plot([0, sigmaNMohr], [C0, tauMohr], '-r', 'LineWidth', 1)  
#hold on ; 
fill([0, Pf, Pf, 0], [0, 0, sigma1, sigma1], 'b', 'FaceAlpha', 0.2) 
contourf(sigmaN, tau, muOAfracture, ncontours, 'EdgeColor', 'none') 
plotMohr_OA(sigma1, sigma2, sigma3, 1, '-k')  
for f in n_fractures: #= 1:n_fractures
    plot(sigmaNFracture(f), tauFracture(f), '.r', 'MarkerSize', sizePoleMarker, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')  
#hold off ; 
colormap(cbOA) 
xlim([0, sigma1*1.05])  
ylim([0, sigmad*0.75]) 
xlabel('Effective normal stress, MPa') 
ylabel('Shear stress, MPa') 
title('Opening angle \mu_a for P_f=', num2str(Pf), ' MPa  n=', num2str(n_fractures)) 
caxis([-10, 90])
cb = colorbar  
cb.Location = 'SouthOutside'  
cb.Label.String = 'Opening angle, \circ'  
#print -r600 -dtiff 'FracTend_OA_mohr.tif'  
 
print(' ') 
print(['*** ...finished FracTend version ', num2str(version), ' at ', datestr(now), '.']) 
print(' ') 
