function [svc_final] = process_svc_3Coutput(input_struct,input_path,station_ids)
% %%%%%%%%%%%%%%%%%%%%%% %
% process_svc_3Coutput.m %
% %%%%%%%%%%%%%%%%%%%%%% %
%
% Inputs:
% -------
%          input_struct : structure containing svc data with 
%                         scans sorted correctly and outliers identified
%                         (default = svc_final [output from 
%                         process_svc_rrs.m])
%           input_path : location of files generated with 3C (in Python)
%          station_ids : Station name/ID for colocated data, aligned with
%                        the input structure
%
% Outputs:
% --------
%             svc_final : structure array from process_svc_rrs, updated
%                         with Rrs values calculated using 3C
%
% -------------------------------------------------------------------------

% Set inputs to default values if not defined
% -------------------------------------------

if exist('input_struct','var')==0
    disp('Input structure required to run function');
    return
end

if exist('input_path','var')==0
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

% Get 3C output file names from input directory
% ---------------------------------------------
contents = dir(input_path);
files = char([]);

for ii = 1:length(contents)

    if ~strcmp(contents(ii).name(1),'.') % removes "hidden" contents
        files = [files; char(contents(ii).name)];
    end

end

% Get station ids from file names (everything before the last underscore)
% -----------------------------------------------------------------------
file_ids = string([]); 

for ii = 1:size(files,1)

    uidx = strfind(files(ii,:),'_');
    file_ids(ii) = string(files(ii,1:max(uidx)-1));

end

file_ids = file_ids';
files = string(files);

clear contents uidx;

% Initiate output structure & append data within 3C output files
% --------------------------------------------------------------

fn = fieldnames(input_struct);

for ii=1:length(input_struct)

    svc_final(ii).station_id = station_ids(ii);

    for jj = 1:length(fn)
        svc_final(ii).(fn{jj}) = input_struct(ii).(fn{jj});
    end

    idx = find(strcmp(file_ids,svc_final(ii).station_id));

    for jj = 1:length(idx)

        dat=importdata(fullfile(input_path,files{idx(jj)}),',');
        svc_final(ii).rrs_wave_3C(:,jj)=dat.data(:,1);
        svc_final(ii).rrs_3C(:,jj)=dat.data(:,5);

    end
    
end

end