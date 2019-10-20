# process_svc
Function to process raw Spectra Vista Corporation radiance into remotely-sensed reflectance.

The **process_svc** function reads in raw above-water radiometry data, sorts spectra into a unique station based on a time threshold (e.g. 15 minutes), clusters spectra to automatically classify as reference, sky or water spectra, identifies outliers and discards, then calculates mean radiance for each scan type. A flag is generated to indicate suspect spectra that do not meet an expected number of quality spectra (e.g., 3 out of 5 spectra identified as good). Using QAQC’d data, remotely-sensed reflectance is calculated using either a user defined or default rho value from *Mobley et al. 1999*. Remotely-sensed reflectance with an applied offset is also calculated based on a user defined or default wavelength range, where remotely-sensed reflectance is assumed to be 0 within this wavelength range, and the average remotely-sensed reflectance value is estimated for this wavelength range and subtracted from the entire spectra. Final output is a Matlab data structure.

Currently, the function auto-populations up to 15 columns of data. If the data file has more columns of data than this, and the wavelength or radiance data is in a column number greater than 15, the user will need to add additional variable names to the ‘col_names’ variable within the function. To add names, users will simply add “ ,’Var#’ “ to the character string, where # indicates the number. So, the first number added would be ‘Var16’, as in the example below:

    col_names={'Var1','Var2','Var3','Var4','Var5','Var6','Var7','Var8','Var9','Var10','Var11','Var12','Var13','Var14','Var15'};

becomes

    col_names={'Var1','Var2','Var3','Var4','Var5','Var6','Var7','Var8','Var9','Var10','Var11','Var12','Var13','Var14','Var15',’Var16};

Users will be notified of an issue with the error message printed to the Matlab console, lines 91-97 of the **process_svc** function.

## Input and Output

**Input**

data_folder          = path to data folder  containing raw svc scan files (required input)
wavelength_column    = column number for wavelength values (e.g., 1)
radiance_column      = column number for radiance values (e.g., 5) 
time_diff            = number of minutes used to define a station
qf_lower             = lower threshold used to flag potential quality issues if fewer scans for a particular scan type are acquired      
qf_upper             = upper threshold used to flag potential quality issues if more scans for a particular scan type are acquired 
rho                  = constant accounting for skylight (0.028=standard, Mobley 1999)
Rg                   = plaque reflectance value (0.99)
min_null_lam         = minimum wavelength for range used to null correct reflectance (e.g., 900 nm in coastal waters)
max_null_lam         = maximum wavelength for range used to null correct reflectance (e.g., 950 nm in coastal waters)

**Output**

*svc_output* = structure containing fields as follows:

latitude             = average latitude across all scans defined for a single station (from time threshold)
longitude            = average longitude across all scans defined for a single station (from time threshold)
quality_flag         = 0 if scans passed quality test, 1 if scans failed quality test
water_lam            = water scan wavelengths
water_rad            = water scan radiance values
water_files          = water source file name(s)
sky_lam              = sky scan wavelengths
sky_rad              = sky scan radiance values
sky_files            = sky source file name(s)
ref_lam              = ref scan wavelengths
ref_rad              = ref scan radiance values
ref_files            = ref source file name(s)
Rrs_lam              = wavelengths corresponding to calculated Rrs values
Rrs                  = Rrs calculated using Eq. 6, Mobley 1999
Rrs_offset           = Rrs calculated using Eq. 7, Mobley 1999

## Example Code to Run Function

    data_folder='~/MyComputer/MyData/SVC_data/scans/';
 
    svc=process_svc(data_folder,15,3,5,0.028,0.99,900,950);

## Example Code to Plot Results

    figure;
    plot(svc(1).Rrs_lam,svc(1).Rrs,'-k')
    set(gca,'XLim',[min(svc(1).Rrs_lam) max(svc(1).Rrs_lam)],'fontsize',16,'fontname','times new roman')
    xlabel('Wavelength (nm)')
    ylabel('R_{rs} (sr^{-1})')

## Example Code to Select All Spectra from July

Note: svc is a structure, so above, we’re calling the first element of that structure to plot. We can make this more complicated by:

    % find all data from July and plot
    for ii=1:length(svc)
        if month(svc(ii).datetime)==7 && svc(ii).quality_flag==0
            ind(ii)=1;
        else
            ind(ii)=0;
        end
    end

    ind=find(ind==1);

    figure;
    for ii=1:length(ind)
        hold on;
        plot(svc(ind(ii)).Rrs_lam,svc(ind(ii)).Rrs)
    end
    set(gca,'XLim',[min(svc(ind(1)).Rrs_lam) max(svc(ind(1)).Rrs_lam)],'fontsize',16,'fontname','times new roman')
    xlabel('Wavelength (nm)')
    ylabel('R_{rs} (sr^{-1})')
