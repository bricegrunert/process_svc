function process_svc_3Cinput(input_struct,output_path,station_ids,Rg,cruise,region,secchi,wind_spd,wind_dir,depth)
% %%%%%%%%%%%%%%%%% %
% process_svc_3Cinput.m %
% %%%%%%%%%%%%%%%%% %
%
% Inputs:
% -------
%          input_struct : structure containing svc data with 
%                         scans sorted correctly and outliers identified
%                         (default = svc_final [output from 
%                         process_svc_rrs.m])
%    %%%               dat : table template for output files
%           output_path : location to put files used as 3C input
%           station_ids : station name for each station with corresponding
%                         radiometry data
%                    Rg : plaque reflectance value (default = 0.99 for 
%                         Spectralon)
%                cruise : name of cruise observations were a part of,
%                         provided as a character array
%                region : general location of observations, provided as a 
%                         character array
%                secchi : observed secchi depth
%              wind_spd : observed wind speed
%              wind_dir : observed wind direction
%                 depth : observed depth
%
% Outputs:
% --------
% 
% No outputs are generated within Matlab. Outputs are .csv files generated
% as input for the three-component (3C) model of Groetsch et al. (2017)
%
% -------------------------------------------------------------------------


% Set inputs to default values if not defined
% -------------------------------------------

if exist('input_struct','var')==0
    disp('Input structure required to run function');
    return
end

if exist('output_path','var')==0
    disp('No path provided for 3C input files');
    return
end


% Number of stations
if exist('station_ids','var')==0
    cnt = 0;
    for ii = 1:length(input_struct)
        cnt = cnt + 1;
        station_ids{ii} = strcat(['St',num2str(cnt)]);
    end
end

% Rg
if exist('Rg','var')==0
    Rg = 0.99;
end

% Cruise
if exist('cruise','var')==0
    cruise = repmat({'Not provided'},length(input_struct),1);
end

% Region
if exist('region','var')==0
    region = repmat({'Not provided'},length(input_struct),1);
end

% Secchi
if exist('secchi','var')==0
    secchi = repmat({'NaN'},length(input_struct),1);
else
    secchi = num2str(secchi);
end

% Air temperature
if exist('at','var')==0
    at = repmat({'NaN'},length(input_struct),1);
else
    at = num2str(at);end

% Sea surface temperature
if exist('sst','var')==0
    sst = repmat({'NaN'},length(input_struct),1);
else
    sst = num2str(sst);
end

% Wind speed
if exist('wind_spd','var')==0
    wind_spd = repmat({'NaN'},length(input_struct),1);
else
    wind_spd = num2str(wind_spd);
end

% Wind direction
if exist('wind_dir','var')==0
    wind_dir = repmat({'NaN'},length(input_struct),1);
else
    wind_dir = num2str(wind_dir);
end

% Depth
if exist('depth','var')==0
    depth = repmat({'NaN'},length(input_struct),1);
else
    depth = num2str(depth);
end


% (6) Create files for input into 3C
% ----------------------------------

wave_3C=350:1:900;

for ii=1:length(input_struct)

    dat{1,1}=['#'];
    dat{2,1}=['#'];
    dat{3,1}=strjoin(['# ID:',station_ids(ii)],' ');
    dat{4,1}=strjoin(['# Mission:',cruise(ii)],' ');
    dat{5,1}=strjoin(['# Location: ' region(ii)],' ');
    dat{6,1}=['# Latitude: ' num2str(input_struct(ii).latitude)];
    dat{7,1}=['# Longitude: ' num2str(input_struct(ii).longitude)];
    dat{8,1}=['# DateTime: ' datestr(input_struct(ii).datetime,'mm/dd/yyyy HH:MM:SS')];
    dat{9,1}=strjoin(['# Secchi Depth [m]: ' secchi(ii)],' ');
    dat{10,1}=strjoin(['# Air Temperature [degC]: ' at(ii)],' ');
    dat{11,1}=strjoin(['# Sea Temperature [degC]: ' sst(ii)],' ');
    dat{12,1}=strjoin(['# Wind Speed [m/s]: ' wind_spd(ii)],' ');
    dat{13,1}=strjoin(['# Wind Direction [deg]: ' wind_dir(ii)],' ');
    dat{14,1}=strjoin(['# Depth [m]: ' depth(ii)],' ');
    dat{15,1}=['#'];    
    dat{16,1}=['Wavelength, [nm]'];    
    dat{16,2}=['Sky Radiance, [W/(m^2 nm sr)]'];    
    dat{16,3}=['Upwelling Radiance, [W/(m^2 nm sr)]'];    
    dat{16,4}=['Downwelling Irradiance, [W/(m^2 nm)]'];    

    mn=nanmean(input_struct(ii).sky_rad,2);
    sky=interp1(input_struct(ii).sky_wave(:,1),mn,wave_3C,'linear');

    mn=nanmean(input_struct(ii).ref_rad,2);
    reference=interp1(input_struct(ii).ref_wave(:,1),mn,wave_3C,'linear');
    irradiance=pi.*reference*(1/Rg);

    cnt = 0;
    for jj = 1:size(input_struct(ii).wat_wave,2)

        cnt = cnt + 1;
        water = interp1(input_struct(ii).wat_wave(:,1),input_struct(ii).wat_rad(:,jj),wave_3C,'linear')';

        dat(17:length(sky)+16,1)=num2cell(wave_3C);
        dat(17:length(sky)+16,2)=num2cell(sky);
        dat(17:size(water,1)+16,3)=num2cell(water);
        dat(17:length(irradiance)+16,4)=num2cell(irradiance);
        clear T filename
        filename=strcat(station_ids(ii),'_',num2str(cnt),'.csv');
        T = array2table(dat);
        writetable(T,fullfile(output_path,char(filename)),'WriteVariableNames',0);

    end

end

end