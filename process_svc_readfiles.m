function [svc_data] = process_svc_readfiles(input_path,num_stations)
% %%%%%%%%%%%%%%%%% %
% process_svc_readfiles.m %
% %%%%%%%%%%%%%%%%% %
%
% Inputs:
% -------
%          input_path : path to raw SVC files (converted to .txt from .sig
%                       format) 
%        num_stations : number of stations where radiometry observations
%                       were recorded
%
% Outputs:
% --------
%            svc_data : structure array containing location/time info and 
%                       plaque, sky, and water radiance data for each 
%                       unique station. Scan types are automatically 
%                       classified using k-means clustering (manually check
%                       classifications for accuracy)
%
% -------------------------------------------------------------------------

% Set inputs to default values if not defined
% -------------------------------------------

if exist('input_path','var')==0
    disp('Input path required to run function');
    return
end

if exist('num_stations','var')==0
    disp('Number of stations with observations is needed');
    return
end


% Get raw file names from input directory
% ---------------------------------------
contents = dir(input_path);
files = string([]);

for ii = 1:length(contents)

    if ~strcmp(contents(ii).name(1),'.') % removes "hidden" contents
        files = [files; string(contents(ii).name)];
    end

end

clear contents;

% Extract datetime, latitude, longitude, wavelength, and radiance data 
% --------------------------------------------------------------------
for ii = 1:length(files)

    % Find number of header lines - provide an "overly generous" starting
    % value, then search for "data="

    tst = importdata(fullfile(input_path,files(ii)),' ',100);

    nhdr = find(contains(tst.textdata(:,1),'data=')==1);

    data = importdata(fullfile(input_path,files(ii)),' ',nhdr);
    dts = string(split(data.textdata(18))); % assumes time info is on row 18
    dtime(ii) = datetime(strcat(dts(5)," ",dts(6)," ",dts(7)),'InputFormat','MM/dd/yyyy hh:mm:ss a');
    lats = char(split(data.textdata(20))); % assumes latitude info is on row 20

    % location of latitude can change, but the one you want is always after
    % the comma
    ind = find(strcmp(split(data.textdata(20)),',')==1) + 1;

    if size(lats,1) < 4
        latitude(ii) = NaN;
    else
        latitude(ii) = str2double(lats(ind,1:2)) + str2double(lats(ind,3:9))/60; % converts from decimal minutes to decimal degrees
        if strcmp(lats(ind,end),'S')
            latitude(ii) = -latitude(ii); % Makes southern latitudes negative
        end
    end

    lons = char(split(data.textdata(19))); % assumes longitude info is on row 19

    % location of latitude can change, but the one you want is always after
    % the comma
    ind = find(strcmp(split(data.textdata(19)),',')==1) + 1;

    if size(lons,1) < 4
        longitude(ii) = NaN;
    else
        longitude(ii) = str2double(lons(ind,1:3)) + str2double(lons(ind,4:10))/60;
        if strcmp(lons(ind,end),'W')
            longitude(ii) = -longitude(ii); % Makes western longitudes negative
        end
    end

    wavelength(:,ii) = data.data(:,1);
    radiance(:,ii) = data.data(:,3);  

end

clear data dts lats lons;


% Identify stations based on date/time & location clusters
% ------------------------------------------------------

clear idx
idx(:,1) = datenum(dtime);
idx(:,2) = latitude;
idx(:,3) = longitude;

[icluster, centroids] = kmeans(idx,num_stations);

[B,I] = sort(centroids(:,1),1);

for ii = 1:length(I)

    ind = find(icluster == I(ii));

    svc_data(ii).datetime = median(dtime(ind),2);
    svc_data(ii).dt_all = dtime(ind);
    svc_data(ii).lat_all = latitude(ind);
    svc_data(ii).latitude = nanmedian(latitude(ind),2);
    svc_data(ii).lon_all = longitude(ind);
    svc_data(ii).longitude = nanmedian(longitude(ind),2);
    svc_data(ii).wavelength = wavelength(:,ind);
    svc_data(ii).radiance = radiance(:,ind);
    svc_data(ii).files = files(ind);

end


% Separate reference, sky, and water scans using k-means clustering
% -----------------------------------------------------------------

for ii = 1:length(svc_data)

    Q = trapz(svc_data(ii).radiance,1);

    istart(1) = min(Q);
    istart(2) = median(Q);
    istart(3) = max(Q);
    [cid,C] = kmeans(Q',3,'start',istart');

    ref_idx = find(cid == find(C == max(C)));
    sky_idx = find(cid == find(C == median(C)));
    wat_idx = find(cid == find(C == min(C)));

    svc_data(ii).ref_files = svc_data(ii).files(ref_idx)';
    svc_data(ii).ref_wave = svc_data(ii).wavelength(:,ref_idx);
    svc_data(ii).ref_rad = svc_data(ii).radiance(:,ref_idx);
    svc_data(ii).sky_files = svc_data(ii).files(sky_idx)';
    svc_data(ii).sky_wave = svc_data(ii).wavelength(:,sky_idx);
    svc_data(ii).sky_rad = svc_data(ii).radiance(:,sky_idx);
    svc_data(ii).wat_files = svc_data(ii).files(wat_idx)';
    svc_data(ii).wat_wave = svc_data(ii).wavelength(:,wat_idx);
    svc_data(ii).wat_rad = svc_data(ii).radiance(:,wat_idx);

end

end
