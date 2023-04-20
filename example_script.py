import numpy as np
import matplotlib.pyplot as plt
from overlay_spectral_lines import SpectralLineOverlays

'''
Author: A. J. Hedge
Created: 5/04/2023
Last Modified: 19/04/2023
Usage: python3 example_script.py
Purpose: Read example spectrum file, demonstrate with a few redshifts the SpectralLineOverlays' use of the catalogue file of
        molecular vibrational emission lines to indicate where they would appear in the spectrum at each redshift solution.
'''

# Setup some parameters (most will just be the class default)
catalogue_file = "vibrational_line_catalogue.dat"
annotate=True
name = 'Test spectra'

# Load example spectrum data and (optionally) plot onto an axes that can be given to the Overlayer
ax = plt.axes()
data = np.genfromtxt("example_spectra.spec")
ax.plot(data[:,0],data[:,1],lw=0.6)
ax.set_xlim(np.min(data[:,0]),np.max(data[:,0]))

# Let's try a few z values
z_data = [0, 1.5, 3, 4.5, 6]

# Initialise the Overlayer
slo = SpectralLineOverlays(catalogue_file,ax=ax,annotate=annotate)

for z_in in z_data:
    ##### DEMO FOR RE-USING AXES USED TO INITIALISE OVERLAYER #####
    ax1 = slo.overlay_all_lines(z=z_in,title_name=name)
    plt.savefig("spectra/"+name+f".z_{z_in:.3f}_lines_keep_ax_demo.png",dpi=100)

    ##### DEMO FOR PASSING IN FREQ & FLUX DATA WITH NO AXES (uses in-built plot_spectrum) #####
    ax2 = slo.overlay_all_lines(z=z_in,freq_data=data[:,0],flux_data=data[:,1],title_name=name)
    plt.savefig("spectra/"+name+f".z_{z_in:.3f}_lines_ax_from_data_demo.png",dpi=100)

    ##### DEMO FOR PASSING IN A NEW AXES EVERY TIME #####
    new_ax = plt.axes()
    new_ax.plot(data[:,0],data[:,1],lw=0.6)
    new_ax.set_xlim(np.min(data[:,0]),np.max(data[:,0]))
    ax3 = slo.overlay_all_lines(new_axes=new_ax,z=z_in,title_name=name)
    plt.savefig("spectra/"+name+f".z_{z_in:.3f}_lines_new_ax_demo.png",dpi=100)

