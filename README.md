# Description
Contains a simple class `SpectralLineOverlays` that overlays redshifted spectral lines onto your own spectra. This code is designed to be imported into your projects and flexibly integrate with your own plotting code by receiving your plot's axes and overlaying rotational emission lines of molecules onto the axes (assuming the molecules behave as rigid rotors). Optionally, individual emission lines specified in an extra file can also be added to the overlay.

# Background
Some key equations that you should understand before using the code are to do with:
- redshifting of emitted light
  * given a redshift $z$, the observed frequency $f_{obs}$ of light at an emitted frequency $f$ is given as:
  
```math
f_{obs}=\frac{f}{1+z}
```
                        
- transition frequency (energy) of adjacent energy levels of a **rotating linear molecule as a rigid rotor**
  * one way of writing the energy of a rigid rotor with angular momentum quantum number $J$ is:
  
```math
E_J=BJ(J+1)
```
    
   where $B$ is the rotational constant of the molecule. $B$ expressed in wavenumbers is $\frac{h}{8\pi^2 cI}$, where I
   is the moment of inertia of the molecule. The moment of inertia is related to the bond length and reduced mass of a
   molecule as $I=\mu l^2$ and remains constant for a single molecule, hence $B$ will also remain constant (for a rigid rotor).
   As such we can write the energy emitted as $hf$ and express the frequency of light emitted by a rigid rotor
   from a $J\to J-1$ transition,
    
$$\begin{align}
E_{J\to J-1}&=E_J-E_{J-1} \\\\\\
hf&=\left(BJ(J+1)\right)-\left(B(J-1)J\right) \\\\\\
&=B(2J) \\\\\\
f_{J\to J-1}&=\frac{B}{h}(2J)
\end{align}$$
                    
   So adjacent energy level transitions of a rigid rotor simply have an emission frequency proportional to $2J$. i.e.
    
$$f_{J\to J-1}\propto 2J$$
                    
   This is useful for sampling a range of redshifts using a single rigid rotor molecule as its emission lines will appear
   at increasing, predictable intervals,
    
$$\frac{f_{J+1\to J}}{f_{J\to J-1}}=\frac{J+1}{J}$$
                    
   and we can easily build up these emission frequencies from the $J=1\to J=0$ frequency (denoted by $f_0$). e.g.
    
$$\frac{f_{J+1\to J}}{f_0}=\frac{J+1}{1}$$
                    
   Hence the emission frequency of any $J+1\to J$ transition relative to the first transition frequency is,
    
$$f_{J+1\to J}=(J+1)f_0$$


# File Structure
The required file structure to use this code is defined by the relative path to the product directory that you define.

`overlay_spectral_lines.py` is the "starting point" and you define the: `relative/path/to/products/` which will contain,
- `spectra/` directory for output figures by the code (to be implemented, save in your own script for now)
- `relative/path/to/line_catalogue` files (the string you provide as the filename can be just the filename if it is in the
  products directory, or a relative path *from* the products directory followed by the filename.)
  * the contents of your `line_catalogue` files must contain columns of data formatted as: \<annotation name> \<frequency> \<matplotlib colour>
  so that it is parsed correctly by the Overlayer (see example `rotational_line_catalogue.dat` which works with the example code).

It is expected that you will read in your own spectrum data wherever it may be in your preferred way, and either format
it to be compatible when providing it to the Overlayer, or create your own base spectrum plot to provide to the Overlayer.

# Usage
You may want to try running the `example_script.py` file and understand how it is interfacing with the Overlayer instance
in each of the demos. A short explanation is as follows,

Import as: `from overlay_spectral_lines import SpectralLineOverlays`

You can then initialise an instance of the Overlayer using the class definition (supplying as many parameters as necessary
for your needs): `slo = SpectralLineOverlays(params...)`

You can then use the following accessible methods from the Overlayer:
- `fn_transitions`: returns the redshifted rotational emission lines in the desired frequency range
- `fs_transitions`: returns the redshifted specified frequency
- `plot_wrapper`: specify the appropriate `f?_transitions` function for a molecule in the lines catalogues (after slicing by
                n_lines) and it will return an axes with the line(s) overlaid (you may optionally provide spectrum data
                OR an already set-up axes with the spectrum)
- `overlay_all_lines`: for a specified redshift, will overlay all the molecular emission lines in the lines catalogues
                     (after slicing by n_lines) and return an axes (you may optionally provide spectrum data OR an
                     already set-up axes with the spectrum)

You can also edit `main` in `overlay_spectral_lines.py` if you instead wish to run from the file with the class definition
directly instead of importing. The usage in this case should be similar to that of importing and then using the class.

# Additions to Future Versions
Here I detail planned functional additions to the code. Bugfixes will not be included here and any "bugs" should be reported
on the [Issues](https://github.com/aj-hedge/overlay_spectral_lines/issues) page. I invite suggestions there as well, but
cannot guarantee I will get round to them!

- [ ] Changeable unit scale (hard-coded GHz and mJy right now)
- [ ] show_figures, save_figures implementation
- [ ] better annotation height algorithm
- [x] More flexible handling of line definitions read from the line_catalogue_file and extra_lines_file
     (e.g. if the user wants only extra_lines_file)

# Acknowledging this work
As per the BSD 3-Clause License, you may not use my name to *endorse or promote* products derived of this work, but I do
request that you provide the following acknowledgement:

This work made use of the code found at https://github.com/aj-hedge/overlay_spectral_lines (Hedge, A. J.) for producing
line overlays on plots.

In future if this original code is included in an accepted paper, I will update this section to include a reference to
the paper that I will also request be cited for derivative works beyond that date.

# References
Sample frequency transitions were obtained from the [ALMA Splatalogue](https://www.splatalogue.online/) with reference to the following line list data,

CDMS: H. S. P. Müller, F. Schlöder, J. Stutzki, and G. Winnewisser, [J. Mol. Struct. 742, 215-227 (2005)](http://dx.doi.org/10.1016/j.molstruc.2005.01.027)

JPL: H. M. Pickett, R. L. Poynter, E. A. Cohen, M. L. Delitsky, J. C. Pearson, and H. S. P. Muller, "Submillimeter, Millimeter, and Microwave Spectral Line Catalog," J. Quant. Spectrosc. & Rad. Transfer 60, 883-890 (1998).
