#AuNP-UVVis-Plot
This is a python script that I used to plot UV-Vis spectra of gold nanoparticles (AuNPs) taken by a *Cary 60* UV-Vis spectrometer (Agilent). The script was not intended for publication, therefore the code is rather messy and only sparingly documented. Maybe it comes handy anyway.

It is likely that this will break if you use CSV files by different devices, so you might have to fix the import logic. Feel free to send pull-requests or to message me if you have questions :)

### Features

- Plot one or multiple UV-Vis spectra
- Command line interface to select CSV files and (if applicable) one of multiple samples within a CSV file.
- Normalize spectra.
- Find and fit the main LSPR peak at ~520 nm wavelength.
- Calculate the concentration of AuNPs based on a given molar extinction coefficent. (You can find them tabulated by [Haiss *et al.*](https://doi.org/10.1021/ac0702084 "Haiss *et al.*"), in the SI)

### Documentation
Basic documentation will be added soon.
