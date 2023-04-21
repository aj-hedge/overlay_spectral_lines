import numpy as np
import matplotlib.pyplot as plt
import pickle
import itertools
from typing import Callable
from types import NoneType
import sys,os,glob

'''
Author: A. J. Hedge
Created: 5/04/2023
Last Modified: 20/04/2023
Usage: Import into script and use non-static functions with own spectrum plotting and line emission catalogue file, or edit main() and run.
Purpose: Take in spectrum files, get the associated redshift-solutions files and use knowledge bank of (assuming simple harmonic oscillator)
        molecular vibrational emission lines to indicate where they would appear in the spectrum at each redshift solution.
'''

# Frequency transitions obtained from the ALMA Splatalogue and are in GHz (cut-off CDMS/JPL intensity > -7)

# TODO:
#       - Changeable unit scale (hard-coded GHz and mJy right now)
#       - show_figures, save_figures implementation
#       - better annotation height algorithm?
#	    - More flexible handling of line definitions read from the line_catalogue_file and extra_lines_file
#		(e.g. if the user wants only extra_lines_file


class SpectralLineOverlays(object):
    def __init__(self, line_catalogue_file: str, ax: plt.Axes=None, extra_lines_file: str=None, product_directory: str='./', n_lines: int=2,
                 annotate: bool=False, spectra_only: bool=False, show_figures: bool=False, save_figures: bool=False, force_quiet: bool=False):
        '''
        `LineOverlays` receives (at least) the name of a file containing a list of molecules' fundamental vibrational transition frequencies
        which it will then use to overlay redshifted counterparts to those frequencies on an axes. Although some parameters have default
        values, it is strongly recommended to change the following parameters to suit your needs/workspace setup.

        Parameters
        ----------
        line_catalogue_file (required) : `str`
            Name of the file containing (in order): molecule_name first_vibrational_transition_frequency matplotlib_colour_string
            columns of data. The contents of this file will be parsed expecting the above format and define the appearance of overlayed
            lines.
        
        ax : `matplotlib.pyplot.Axes`
            The axes instance that the overlays will be placed on. If left as `None`, an included function will be called to plot a
            random spectra and provide the axes instance. Providing a new axes or freq & flux data later to plotting wrapper functions
            will replace the placeholder axes and likely resolve any issues that would otherwise arise. If you are iteratively
            making plots with overlays, the overlay_all_lines function will clean-up lines & annotations for you!
        
        extra_lines_file : `str`
            Name of the file containing (in order): molecule_name specific_transition_frequency matplotlib_colour_string.
            The contents of this file will be parsed expecting the above format and define the appearance of overlayed lines. This extra
            file is for inserting specific molecules/frequency transitions that DON'T repeat as n*first_vibrational_transition_frequency.
        
        product_directory : `str`
            This is the relative path (from this python file) to the directory containing any required input files and the spectra output
            directory. If the spectra directory is not found at the given product_directory (default is at this file's directory), then
            the spectra directory will be created automatically.
        
        n_lines : `int`
            If you wish to limit (or unlimit) the number of lines included in the overlay that have been read from the line_catalogue_file,
            you can do so by setting n_lines. The default is 2, and by choosing `None` all lines will be selected.
        
        annotate : `bool`
            Set to `True` if you want molecule names annotated next to their respective lines on the overlay (can get cluttered!).
        
        spectra_only : `bool`
            Similar to when ax=None, this option when set to `True` uses the ungeneralised spectrum plotting function and does not add
            line overlays to it. Not intended for use when the user has provided an axes instance.
        
        show_figures : `bool`
            Enables matplotlib.pyplot.show() function call, displaying figures sequentially.
        
        save_figures : `bool`
            Automatically saves the figures once the lines have been overlayed to the spectra directory. Leave as false if you want to
            make manual adjustments to the returned axes before saving it yourself.
        '''

        # Main data
        self.line_catalogue_file = line_catalogue_file
        self.ax = ax
        self.init_ax = pickle.dumps(self.ax)  # Necessary for deep-copy backup
        self.extra_lines_file = extra_lines_file
        self.product_directory = product_directory
        self.n_lines = n_lines
        self.annotate = annotate
        self.spectra_only = spectra_only
        self.show_figures = show_figures
        self.save_figures = save_figures
        self.quiet = force_quiet
        # Tweakable parameters
        self.vline_alpha = 0.6
        self.vline_lw = 1
        self.cat_header_lines = 0

        # Initialise N-D data structures
        self.vibrational_molecules = {}
        self.extra_molecules = {}

        # Navigating and setting up directory
        if os.path.isdir(self.product_directory):
            os.chdir(self.product_directory)
        else:
            print(f"No such relative work directory '{self.product_directory}'")
            exit()

        if os.path.isdir('spectra') == False:
            print("Output directory 'spectra' not found, creating it now.")
            os.mkdir('spectra')

        # Read files and populate dictionaries
        self.vibrational_molecules = self.read_line_catalogue(self.vibrational_molecules, self.line_catalogue_file,self.cat_header_lines)
        if self.extra_lines_file != None:
            self.extra_molecules = self.read_line_catalogue(self.extra_molecules, self.extra_lines_file,self.cat_header_lines)
        
        # Avoid unnecessary index out of bounds errors when the user clearly wants all lines (but specified more)
        if n_lines == None:
            n_lines = len(self.vibrational_molecules.keys())
        if n_lines > len(self.vibrational_molecules.keys()):
            n_lines = len(self.vibrational_molecules.keys())
        # if self.extra_lines_file != None and n_lines * 2 > self.extra_molecules.keys():
        #     n_lines = np.floor(len(self.extra_molecules.keys())/2)

        # Slice dictionaries by n_lines if not None
        self.vibrational_molecules = dict(itertools.islice(self.vibrational_molecules.items(),0,n_lines))
        # if self.extra_lines_file != None:
        #     self.extra_molecules = dict(itertools.islice(self.extra_molecules.items(),0,2*n_lines))

        # If no axes given, create a randomised dummy spectrum
        if self.ax == None:
            self.ax = plt.axes()
            freq = np.arange(50,150,0.1)
            flux = np.random.randn(len(freq))
            source_name = 'Placeholder'
            self.ax, self.init_ax = self.plot_spectrum(self.ax,freq,flux,source_name)


    # Static methods
    @staticmethod
    def blockPrint(quiet):
        if quiet:
            sys.stdout = open(os.devnull,'w')

    @staticmethod
    def enablePrint(quiet):
        if quiet:
            sys.stdout = sys.__stdout__

    @staticmethod
    def read_line_catalogue(out_dict: dict, fname: str, n_headerlines: int=0):
        '''
        Reads given line catalogue file and builds the dictionary in the Overlayer instance.

        Parameters
        ----------
        out_dict : dict
            Dictionary to store the information in.
        fname : str
            Name of line catalogue file.

        Returns
        -------
        out_dict : dict
            Dictionary to store the information in.
        '''
        # Read catalogue file data
        lines = []
        with open(fname,'r') as f:
            lines = f.readlines()
        # Parse catalogue file data
        lines_data = [line.replace('\n','').split(' ') for line in lines[n_headerlines:]]
        for dat in lines_data:
            dat[1] = float(dat[1])
        # Build up dictionary
        transitions = np.array([t[1] for t in lines_data],ndmin=1)
        names = np.array([n[0] for n in lines_data],ndmin=1)
        colours = np.array([c[2] for c in lines_data],ndmin=1)
        for trans, name, col in zip(transitions,names,colours):
            out_dict[name] = [trans, col]

        return out_dict
    
    @staticmethod
    def plot_spectrum(ax: plt.Axes, freq: np.ndarray, flux: np.ndarray, source_name: str=''):
        '''
        Adds basic plot of the spectra to the Overlayer instance's stored axes (strict parameters, units).

        Parameters
        ----------
        ax : matplotlib.pyplot.Axes
            Axes to plot spectrum onto.
        freq : numpy.ndarray (float)
            Frequency vector of spectrum.
        flux : numpy.ndarray (float)
            Flux vector of spectrum.
        source_name : str
            Name to put in plot title.

        Returns
        -------
        ax : matplotlib.pyplot.Axes
            The modified axes from self.
        init_ax_copy : matplotlib.pyplot.Axes
            Deep-copy of ax from pickle dumping and loading (used for reloading initial axes prior to overlaying lines).
        '''
        # Base spectrum plot, replaces current axes. Future version may use ax parameter more flexibly.
        ax.remove()
        ax = plt.axes()
        ax.plot(freq,flux,'k',lw=0.5)
        ax.set_title(source_name)
        ax.set_xlabel("Frequency / GHz")
        ax.set_ylabel("Flux Density / mJy")
        ax.set_xlim([np.min(freq), np.max(freq)])
        init_ax_copy = pickle.dumps(ax)

        return ax, init_ax_copy

    @staticmethod
    def handle_new_axes(curr_axes: plt.Axes, new_axes: plt.Axes, init_axes: plt.Axes, freq_data: np.ndarray, flux_data: np.ndarray,
                        clean_overlay: bool=False, quiet: bool=False):
        '''
        Handles different cases when the user is starting a new plotting action which requires either a new axes, existing axes, or
        a new set of frequency and flux data to pass to the in-built spectrum plotting function. Optional cleaning of existing axes.

        Parameters
        ----------
        curr_axes : matplotlib.pyplot.Axes
            The existing axes in case nothing else is provided.
        new_axes : matplotlib.pyplot.Axes
            The new_axes to use instead.
        init_axes : matplotlib.pyplot.Axes
            Initial axes passed in case the axes need to be reset.
        freq_data : numpy.ndarray
            Vector of frequency data to pass to plot_spectrum if necessary.
        flux_data : numpy.ndarray
            Vector of flux data to pass to plot_spectrum if necessary.
        clean_overlay : bool
            Option to clean the existing axes' overlay (spectral lines and annotations).
        quiet : bool
            Supresses messages, but not Warnings or Errors.
        '''
        if type(new_axes) != NoneType and (type(freq_data) != NoneType and type(flux_data) != NoneType):
            print("WARNING: Both a new axes and freq & flux data provided. Only the new axes will be used.")

        SpectralLineOverlays.blockPrint(quiet)

        if type(new_axes) == NoneType:
            if curr_axes.get_title() == 'Placeholder' and (type(freq_data) == NoneType or type(flux_data) == NoneType):
                raise RuntimeError("No axes provided when initialising Overlayer object and the frequency and/or flux data given " + \
                                   "to plot on an axes instead were either incomplete or missing.")
            elif type(freq_data) != NoneType and type(flux_data) != NoneType:
                if len(freq_data) != len(flux_data):
                    raise ValueError("Length of frequency (X) and flux (Y) data do not match!")
                print("Updating axes with provided frequency and flux data.")
                new_axes, init_axes = SpectralLineOverlays.plot_spectrum(curr_axes,freq_data,flux_data)
            else:
                print("Using existing axes stored in the Overlayer object. If this was not desired, check that the frequency " + \
                      "and/or flux data you provided are not None.")
                if clean_overlay == True:
                    # Clean axes of vlines and annotations
                    curr_axes.remove()
                    new_axes = pickle.loads(init_axes)
                else:
                    new_axes = curr_axes
        else:
            curr_axes.remove()
        
        SpectralLineOverlays.enablePrint(quiet)

        return new_axes

    @staticmethod
    def plot_lines(ax: plt.Axes, f_shifted: float, molecule: str, colour: str,  annotate: bool, n_range: np.ndarray=None, alpha: float=0.6, lw: float=1):
        '''
        Adds the spectral line overlays to the Overlayer instance's stored axes.

        Parameters
        ----------
        ax : matplotlib.pyplot.Axes
            Axes to overlay lines onto.
        f_shifted : float
            The redshifted spectral line frequency from f?_transitions.
        molecule : str
            Name of the molecule associated with the line.
        colour : str
            Name of a matplotlib colour.
        annotate : bool
            Flag for allowing molecule name annotations next to lines.
        n_range : numpy.ndarray (int)
            Range of vibrational energy levels where the emission frequency is close to the spectrum frequency range.
        alpha : float
            Alpha value to pass to axvline.
        lw : float
            Linewidth value to pass to axvline.

        Returns
        -------
        ax : matplotlib.pyplot.Axes
            The modified axes from self.
        '''
        if type(n_range) != NoneType:
            for n, line in zip(n_range,f_shifted):
                ax.axvline(line,alpha=alpha,lw=lw,color=colour)
                if annotate:
                    ax.annotate(molecule+f'({n}-{n-1})',(line,ax.get_ylim()[1]*(line % 10)/10))
        else:
            ax.axvline(f_shifted,alpha=alpha,lw=lw,color=colour)
            if annotate:
                ax.annotate(molecule,(f_shifted,ax.get_ylim()[1]*(line % 10)/10))

        return ax

    # Callable methods
    def fn_transitions(self, f0_transition: float, z: float, f_min: float, f_max: float):
        '''
        Calculates redshifted vibrational frequencies within the given frequency range.

        Parameters
        ----------
        f0_transition : float
            The first vibrational energy level transtion (1->0) emission frequency.
        z : float
            The redshift used to shift frequencies.
        f_min : float
            Lower limit of range (usually from the spectra plotted on the axes)
        f_max : float
            Upper limit of range (usually from the spectra plotted on the axes)
        
        Returns
        -------
        n_range : numpy.ndarray (int)
            Range of vibrational energy levels where the emission frequency is close to the spectrum frequency range.
        fn_shifted : numpy.ndarray (float)
            The redshifted spectral line frequencies corresponding to the n_range transitions.
        '''
        n_upper = int(np.ceil(f_max / f0_transition * (1+z)))
        n_lower = int(np.floor(f_min / f0_transition * (1+z)))
        if n_lower == 0:
            n_lower = 1
        if n_upper < 1:
            raise ValueError
        # to be on the safe side, we'll just pad out the n_range by 1 on either side (not going to zero though)
        if n_lower > 1:
            n_lower -= 1
        n_upper += 1
        n_range = np.arange(n_lower,n_upper+1)
        fn_shifted = n_range * f0_transition / (1+z)

        return n_range, fn_shifted

    def fs_transitions(self, fspec_transition: float, z: float, f_min: float, f_max: float):
        '''
        Calculates redshifted vibrational frequency.

        Parameters
        ----------
        fspec_transition : float
            The emission frequency of a specific transition.
        z : float
            The redshift used to shift frequency.
        f_min : float
            Lower limit of range (usually from the spectra plotted on the axes).
        f_max : float
            Upper limit of range (usually from the spectra plotted on the axes).
        
        Returns
        -------
        None
            In place of n_range so this function can be abstracted with fn_transitions.
        fn_shifted : numpy.ndarray (float)
            The redshifted spectral line frequencies corresponding to the n_range transitions.
        '''
        fspec_shifted = fspec_transition / (1+z)

        return None, fspec_shifted

    def plot_wrapper(self, trans_func: Callable, f_trans: float, z: float, molecule: str, colour: str,
                     new_axes: plt.Axes=None, freq_data: np.ndarray=None, flux_data: np.ndarray=None):
        '''
        Wrapper to call frequency transition function and pass outputs to the plotting function to overlay results.
        Overlays are not removed in this function, so it can be called multiple times to stack lines.

        Parameters
        ----------
        trans_func : function
            Abstracted reference to one of the transition functions (fn_transitions, fs_transitions).
        z : float
            The redshift used to shift frequencies.
        molecule : str
            Name of the molecule associated with the line(s).
        colour : str
            Name of a matplotlib colour.
        new_axes : matplotlib.pyplot.Axes
            Use this parameter to update the axes (with a *new* spectrum plotted on it), otherwise the previously
            saved axes to the Overlayer instance will be used (which is fine if using the same spectrum).
        freq_data : numpy.ndarray (float)
            Optional vector of frequencies if using in-built plot_spectrum function (do not provide new_axes).
        flux_data : numpy.ndarray (float)
            Optional vector of flux data if using in-built plot_spectrum function (do not provide new_axes).
        
        Returns
        -------
        ax : matplotlib.pyplot.Axes
            The modified axes from self.
        '''
        
        # Handle axes provided or generate if required and data is instead provided
        self.ax = self.handle_new_axes(self.ax,new_axes,self.init_ax,freq_data,flux_data,False,self.quiet)
        
        annotate = self.annotate

        x_lims = self.ax.get_xlim()
        n_range, f_shifted = trans_func(f_trans,z,x_lims[0],x_lims[1])
        self.ax = self.plot_lines(self.ax,f_shifted,molecule,colour,annotate,n_range,self.vline_alpha,self.vline_lw)

        return self.ax

    def overlay_all_lines(self, new_axes: plt.Axes=None, z: float=0, freq_data: np.ndarray=None, flux_data: np.ndarray=None, title_name: str=''):
        '''
        Top-level wrapper to iterate through molecule dictionaries and feed parameters into e.g. plot_wrapper.
        Accepts desired redshift solution to use on the axes provided when initialising the instance.
        If using the existing axes in the Overlayer, the overlays will be removed automatically.

        Parameters
        ----------
        new_axes : matplotlib.pyplot.Axes
            Use this parameter to update the axes (with a *new* spectrum plotted on it), otherwise the previously
            saved axes to the Overlayer instance will be used (which is fine if using the same spectrum).
        z : float
            Redshift solution to use when applying line overlays to instance's axes.
        freq_data : numpy.ndarray (float)
            Optional vector of frequencies if using in-built plot_spectrum function (do not provide new_axes).
        flux_data : numpy.ndarray (float)
            Optional vector of flux data if using in-built plot_spectrum function (do not provide new_axes).
        title_name : str
            Label to put in the title followed by z=? solution.
        
        Returns
        -------
        ax : matplotlib.pyplot.Axes
            Axes with finalised overlays.
        '''

        # Handle axes provided or generate if required and data is instead provided
        self.ax = self.handle_new_axes(self.ax,new_axes,self.init_ax,freq_data,flux_data,True,self.quiet)

        self.ax.set_title(title_name+f" z={z:.3f} solution")
        fig = plt.gcf()
        fig.set_size_inches(12,5)

        # overlay lines
        for mol, vals in self.vibrational_molecules.items():
            fn = vals[0]
            colour = vals[1]
            self.ax = self.plot_wrapper(self.fn_transitions,fn,z,mol,colour)
        if self.extra_lines_file != None:
            for mol, vals in self.extra_molecules.items():
                fs = vals[0]
                colour = vals[1]
                self.ax = self.plot_wrapper(self.fs_transitions,fs,z,mol,colour)

        return self.ax


# TEMPLATE MAIN when called as script (fill in appropriate values yourself)
def main():
    line_catalogue_file = "path/from/products/to/catalogue.txt"
    product_directory = "relative/path/to/products/"
    annotate = True
    ax = None

    slo = SpectralLineOverlays(line_catalogue_file,ax=ax,product_directory=product_directory,annotate=annotate)

    filenames = glob.glob("globstring for choosing files containing spectrum data")
    for filename in filenames:
        # REPLACE WITH SOME METHOD OF OBTAINING DATA
        data = np.genfromtxt(filename,skip_header=4)
        # CONVERT VALUES TO GHz and mJy, E.G.
        freq = data[:,2] / 1e3  # convert MHz to GHz
        flux = data[:,4] * 1e3  # convert Jy to mJy

        # REPLACE WITH SOME METHOD OF OBTAINING REDSHIFT SOLUTIONS TO USE
        try:
            z_data = np.genfromtxt(filename+'.redshifts',ndmin=1)   # ndmin=1 required to prevent 0-dimension array raising TypeError in for loop (numpy>=1.23)
        except TypeError as e:
            print(f"TypeError: {e} when reading {filename}.\nNote that ndmin requires numpy>=1.23. A work-around will automatically be applied now.")
            z_data = np.array(z_data,ndmin=1)

        # ITERATIVELY MAKE LINE OVERLAYS AND SAVE
        for z_in in z_data:
            # overlay lines
            slo.ax = slo.overlay_all_lines(z=z_in,freq_data=freq,flux_data=flux,title_name=filename)
            
            # save file
            plt.savefig(product_directory + "spectra/"+filename+f".z_{z_in}_lines.png",dpi=100)


if __name__ == '__main__':
    main()

