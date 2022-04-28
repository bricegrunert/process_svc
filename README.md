# process_svc_v2.0

## Brief Description
This package was built to read in raw above-water radiometry data collected by a Spectra Vista Corporation (SVC) HR-512i spectroradiometer. The package consists of 4 core MATLAB functions and the three component (3C) Python model produced by Groetsch et al. 2017 and available on Philipp Grötsch’s GitLab page (rrs_model_3C). The packages are designed to read in SVC files (.sig) that include raw radiance scans, and through several functions generate quality controlled remote-sensing reflectance (*Rrs*) spectra. Specific components of each function are described in the header information (comments) of each .m file. The four MATLAB functions are as follows:

## process_svc_readfiles.m
Files are read into MATLAB and unique stations are identified through *k*-means clustering of datetime, latitude and longitude data and a user-specified number of stations. Then, radiance scans for each station are sorted using *k*-means clustering of the integral of each scan into reference plaque, sky and water scans. A data structure is generated with each element of the structure corresponding to a unique station.

Users should check that scans are properly sorted, particularly for sub-optimal (cloudy) conditions where sky and reference plaque scans are often of similar magnitude and spectral shape (and thus often intermingled during sorting).

## process_svc_rrs.m
The data structure generated from *process_svc_readfiles.m* is used, along with a rho and Rg value (see Mobley 1999) as input variables. Radiance scans are quality controlled using an outlier analysis, with outliers determined as anything greater than 3 median absolute deviations from the median. Outliers are identified as any spectra with more than 25 wavelengths identified as individual outliers. Quality scans are kept, and ‘Mobley’ Rrs values are calculated using Eq. 6 of Mobley 1999.

## process_svc_3Cinput.m
The data structure generated from *process_svc_rrs.m* is used as input, along with the path (location) to write .csv data files, station IDs and header information used by 3C (e.g., Rg, secchi, wind speed). This function does not provide any data processing or analysis, it only renders input files in the appropriate format for the *3C_process_svc.py* script. An important note on station IDs – this function will run without them, and auto-generate a generic ID in numeric order. We highly recommend providing your own to easily link this data back to your other datasets.

## process_svc_3Coutput.m
The data structure generated from *process_svc_rrs.m* is used as input, along with the location of the output files generated from the *3C_process_svc.py* script and station IDs (again, will resort to generic IDs if none are provided). This function does not provide any data processing or analysis, it only populates the data structure with 3C generated Rrs spectra, matched to the corresponding station ID.

## Inputs and Outputs
Inputs and outputs are defined as part of each function and can be viewed with *help ‘function_name’* in the MATLAB console (be sure to add the functions to your MATLAB directory).

## Notes
The original 3C repository from GitLab is included in the provided files, consistent with its license (GNU Lesser General Public License v3.0). Users will need to set up their Python environment to work with 3C. There are many great resources for how to do this online, and we will not profess to be one of them. We have found Spyder to be the most intuitive for running the Python scripts, consistent with our MATLAB heritage. If you do have questions about basic setup and getting this code working, please do not hesitate to reach out – our contact info is below.

## Questions?
Contact Brice Grunert (b.grunert@csuohio.edu) or Kyle Turner (k.turner00@ccny.cuny.edu).

## Acknowledgements
A huge thank you to Philipp Grötsch for answering too many questions and supporting implementation of this package (and also for being a rad dude).

## References
Groetsch, P. M., Gege, P., Simis, S. G., Eleveld, M. A., & Peters, S. W. (2017). Validation of a spectral correction procedure for sun and sky reflections in above-water reflectance measurements. Optics Express, 25(16), A742-A761.

Groetsch, Philipp. (2017). Python implementation of the 3C Water Surface Reflection Model. Zenodo. https://doi.org/10.5281/zenodo.293851

Mobley, C. D. (1999). Estimation of the remote-sensing reflectance from above-surface measurements. Applied optics, 38(36), 7442-7455.

