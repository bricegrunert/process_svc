%% Example script for process_svc_v2.0 package

% Ensure process_svc-v2.0 is added to your MATLAB environment!

%% process_svc_readfiles
% -------------------------------------------

% Location of the raw spectroradiometer radiance (.sig) files
input_path = '/Users/your_name/Documents/svc_data_folder/raw_data/';

% Number of stations where spectroradiometer data was collected
num_stations = 12;

% Run process_svc_readfiles to assign radiance scans to each station and 
% category (reference plaque, sky, water)
svc_data = process_svc_readfiles(input_path,num_stations);

%% NOTE: 
%       Before proceeding to the next step, you should ensure that scan types 
%       are properly assigned and re-order them manually if necessary.
%       This should only be an issue in sub-optimal (cloudy) conditions.

%% process_svc_rrs
% -------------------------------------------

% Define rho and Rg (from Mobley 1999, Eq. 6)
rho = 0.028;
Rg = 0.99;

% Run process_svc_rrs to quality control radiance scans and calculate
% Mobley Rrs using the 'svc_data' structure as input, creating a new
% structure 'svc_final'
svc_final = process_svc_rrs(svc_data,rho,Rg);

%% process_svc_3Cinput
% -------------------------------------------

% Location to write the 3C input data files to
output_path = '/Users/your_name/Documents/svc_data_folder/3C_input/';

% Run process_svc_3Cinput to create 3C input data files and save them to
% the defined location using the 'svc_final' structure as input
process_svc_3Cinput(svc_final,output_path)

% -------------------------------------------------------------------------
%% Time to run 3C_process_svc.py - be sure to update the paths in that file
% -------------------------------------------------------------------------

%% process_svc_3Cinput
% -------------------------------------------

% Location where the 3C output data files were written to
input_path = '/Users/your_name/Documents/svc_data_folder/3C_output/';

% Run process_svc_3Coutput to add the 3C Rrs to the 'svc_final' structure
process_svc_3Coutput(svc_final,input_path);

