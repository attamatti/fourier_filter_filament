USAGE: fourier_filter_filament.py <mrc file> <pixel size>

Will work best on filaments extracted and straightened with bsoft.  
Filament must be oriented vertically.

Currently hardcoded for actin.
Filters in fourier space to eliminate all signal except a box from 70 to 45 Angstrom in the direction perpendicular to the fibril axis and  infinity to 45 angstrom in the direction parallel to the axis.

Collapses the inverse transform in the x dimension, fits a running average to smooth out the peaks, and finds the crossovers by looking for dips.

Later will make make it find the layerlines on its own and calculate filter parameters on the fly. so it can be used on other filaments.
Although probably won't work on amyloid because there's not enough signal at the 4.8A subunit rise.

20190212 - Matt I
