function [svc_final] = process_svc_rrs(input_struct,rho,Rg)
% %%%%%%%%%%%%%%%%% %
% process_svc_rrs.m %
% %%%%%%%%%%%%%%%%% %
%
% Inputs:
% -------
%          input_struct : structure containing svc data with 
%                         stations and scan types sorted correctly 
%                         (default = svc_data [output from 
%                         process_svc_readfiles.m])
%                   rho : constant accounting for skylight (default = 
%                         0.028, standard from Mobley, 1999)
%                    Rg : plaque reflectance value (default = 0.99 for 
%                         Spectralon)
%
% Outputs:
% --------
%             svc_final : structure array containing plaque, sky, and 
%                         water radiances with outliers removed and the 
%                         remote sensing reflectance calculated using the 
%                         standard Mobley, 1999 method
%
% -------------------------------------------------------------------------

% Set inputs to default values if not defined
% -------------------------------------------

if exist('input_struct','var')==0
    disp('Input structure required to run function');
    return
end

if ~exist('rho','var')
    rho = 0.028;
end

if ~exist('Rg','var')
    Rg = 0.99;
end


% Initiate final svc structure
% ----------------------------

svc_final = input_struct;


% Identify and remove outlier spectra and apply quality flags
% -----------------------------------------------------------

for ii = 1:length(svc_final)

    wave_range_ref = find(svc_final(ii).ref_wave(:,1) >= 400 & svc_final(ii).ref_wave(:,1) <= 600);
    outliers_ref = isoutlier(svc_final(ii).ref_rad(wave_range_ref,:),2);
    noutliers_ref = sum(outliers_ref);
    good_ref = find(noutliers_ref < 25);

    wave_range_sky = find(svc_final(ii).sky_wave(:,1) >= 400 & svc_final(ii).sky_wave(:,1) <= 600);
    outliers_sky = isoutlier(svc_final(ii).sky_rad(wave_range_sky,:),2);
    noutliers_sky = sum(outliers_sky);
    good_sky = find(noutliers_sky < 25);

    wave_range_wat = find(svc_final(ii).wat_wave(:,1) >= 400 & svc_final(ii).wat_wave(:,1) <= 600);
    outliers_wat = isoutlier(svc_final(ii).wat_rad(wave_range_wat,:),2);
    noutliers_wat = sum(outliers_wat);
    good_wat = find(noutliers_wat < 25);

    if length(good_ref) < 3 || length(good_sky) < 3 || length(good_wat) < 3
        svc_final(ii).quality_flag = 1;
    else
        svc_final(ii).quality_flag = 0;
    end

    svc_final(ii).ref_files = svc_final(ii).ref_files(good_ref);
    svc_final(ii).ref_wave = svc_final(ii).ref_wave(:,good_ref);
    svc_final(ii).ref_rad = svc_final(ii).ref_rad(:,good_ref);
    svc_final(ii).sky_files = svc_final(ii).sky_files(good_sky);
    svc_final(ii).sky_wave = svc_final(ii).sky_wave(:,good_sky);
    svc_final(ii).sky_rad = svc_final(ii).sky_rad(:,good_sky);
    svc_final(ii).wat_files = svc_final(ii).wat_files(good_wat);
    svc_final(ii).wat_wave = svc_final(ii).wat_wave(:,good_wat);
    svc_final(ii).wat_rad = svc_final(ii).wat_rad(:,good_wat);

    % Calculate remote sensing reflectance using Mobley (1999) Eq. 6
    % --------------------------------------------------------------

    for jj = 1:size(svc_final(ii).wat_wave,2)

        svc_final(ii).rrs_wave(:,jj) = svc_final(ii).wat_wave(:,jj);
        svc_final(ii).rrs_mobley(:,jj) = (svc_final(ii).wat_rad(:,jj)-(rho.*nanmean(svc_final(ii).sky_rad,2)))./((pi.*nanmean(svc_final(ii).ref_rad,2))*(1/Rg));

    end

end

end
    
    
    
    