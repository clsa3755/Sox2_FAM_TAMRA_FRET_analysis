This code was written for an experiment in which DNA transcription factor, Sox2,
was isolated while bound to DNA on its specific domain.
FRET experimentation was carried out to ultimately determine the bend angle of
DNA with varying lengths of Sox2 protein bound to it.
The donor fluorophore, FAM, and the acceptor fluorophore, TAMRA, were attached
to opposite 5' ends of the dsDNA Sox2 binding region.
The excitation wavelengths for FAM and TAMRA were 490 nm and 560 nm respectively.
This code can be used for analogous experiments with either or both a differing
fluorophore and SOX protein.

Data that was calculated or identified before running through the code is as
follows:
absorbance ratio of εFAM490/εTAMRA560 = 0.543 and is constant
emission ratio of εTAMRA490/εTAMRA560 = 0.121 and is constant
first location of peak donor emission data = 520 nm
second location of peak donor emission data = 580 nm
*If this code is to be used with differing fluorophores or SOX protein, these
values will have to be changed.
That circumstance is reflected in the comments of the code.

The tasks that the code executes are as follows:
-Read donor (FAM) emission spectra
-Read acceptor (TAMRA) emission spectra
-Read background (FAM donor only) intensities at donor emission wavelengths
-Extract all necessary data into workable arrays with floats
-Normalize donor emission spectra
-Calculate Fret Effect, Fret Efficiency, distance between fluorophores
(Rd(a) or r), delta r, and DNA bend angle in radians and degrees
-Plot the following spectra: donor emission spectra, acceptor emission spectra,
fit donor spectra, normalized donor spectra extracted
Plot DNA bend angle vs length of Sox2 protein on complex

To execute the code you will need your data to be in three separate .csv files:
-One file will be the donor (FAM) emission intensities at each wavelength after
excitation for each length of Sox2 protein. Wavelengths will need to be in the
left most column and lengths of Sox2 will have to be the header. This should be
easily formatted from each .csv file obtained from the spectroscopy experiments.
-Another file will have the same formatting but be for the acceptor emission
intensity values. (These wavelengths are different and there are less of them)
-The last file will be the emission intensities from the excitation of the FAM
donor only with all corresponding donor emission wavelengths.

To execute in practice, follow comment directives. All that will need to be changed is the .csv files names if they are different from what the files are named in the example and the calculated constants and peak wavelengths.

Example .csv files that were used to write the example code are included in repository.
