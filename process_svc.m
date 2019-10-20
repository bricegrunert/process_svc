function [svc_output]=process_svc(data_folder,wavelength_column,radiance_column,time_diff,qf_lower,qf_upper,rho,Rg,min_null_lam,max_null_lam)

% daisea
%
% [svc_output]=process_svc(enclosing_folder,data_folder,time_diff,rho,Rg,min_null_lam,max_null_lam)
%
% Function reads in raw Spectra Vista Corporation HR series
% spectroradiometer data, sorts by station and scan type (here, reference, 
% sky and water surfaces) through k means clustering, and discards outlier
% spectra for each type. Remotely-sensed reflectance is calculated using
% the approach of Mobley 1999, with custom or pre-defined constants.
%
% Inputs:
% 
%           data_folder          = path to data folder  containing raw svc 
%                                  scan files (required input)
%           wavelength_column    = column number for wavelength values
%                                  (e.g., 1)
%           radiance_column      = column number for radiance values
%                                  (e.g., 5)
%           time_diff            = number of minutes used to define a 
%                                  station
%           qf_lower             = lower threshold used to flag potential
%                                  quality issues if fewer scans for a 
%                                  particular scan type are acquired
%           qf_upper             = upper threshold used to flag potential
%                                  quality issues if more scans for a 
%                                  particular scan type are acquired
%           rho                  = constant accounting for skylight
%                                  (0.028=standard, Mobley 1999)
%           Rg                   = plaque reflectance value (0.99)
%           min_null_lam         = minimum wavelength for range used to 
%                                  null correct reflectance 
%                                  (e.g., 900 nm in coastal waters)
%           max_null_lam         = maximum wavelength for range used to 
%                                  null correct reflectance 
%                                  (e.g., 950 nm in coastal waters)
%
% Returns:
%     svc_output = structure containing fields as follows:
%
%           latitude             = average latitude across all scans
%                                  defined for a single station 
%                                  (from time threshold)
%           longitude            = average longitude across all scans
%                                  defined for a single station 
%                                  (from time threshold)
%           quality_flag         = 0 if scans passed quality test, 1 if
%                                  scans failed quality test
%           water_lam            = water scan wavelengths
%           water_rad            = water scan radiance values
%           water_files          = water source file name(s)
%           sky_lam              = sky scan wavelengths
%           sky_rad              = sky scan radiance values
%           sky_files            = sky source file name(s)
%           ref_lam              = ref scan wavelengths
%           ref_rad              = ref scan radiance values
%           ref_files            = ref source file name(s)
%           Rrs_lam              = wavelengths corresponding to calculated
%                                  Rrs values
%           Rrs                  = Rrs calculated using Eq. 6, Mobley 1999
%           Rrs_offset           = Rrs calculated using Eq. 7, Mobley 1999
%
%
%
%References
%
%Mobley, Curtis D. "Estimation of the remote-sensing reflectance from 
%above-surface measurements." Applied optics 38.36 (1999): 7442-7455.
%
%
%
% copyright (c) 2019 Brice K. Grunert
% email: bricegrunert@gmail.com
%
% Created 1 February 2019 by BG
%
% Last modified on 17 October 2019 by BG


%%

% code credit AC 2008
d=dir(data_folder);
files={};
[files{1:length(d)}]=deal(d.name);
files=files(setdiff(1:length(files),find(strncmp('.',files,1)))); % get rid of . directories

col_names={'Var1','Var2','Var3','Var4','Var5','Var6','Var7','Var8','Var9','Var10','Var11','Var12','Var13','Var14','Var15'};

if wavelength_column > length(col_names) | radiance_column > length(col_names)
    fprintf('Column number provided exceeds pre-populated column names \n')
    fprintf('User needs to generate corresponding VarX value for column names \n')
    fprintf('Please see ReadMe file for instructions \n')
    svc_output='Invalid Entry';
    return
end


param1=col_names{wavelength_column};
param2=col_names{radiance_column};

for ii=1:length(files)
    dat=readtable(fullfile(data_folder,files{ii}));
    wavelength(:,ii)=dat.(param1);
    radiance(:,ii)=dat.(param2);
end


%% fill in variables if no input

if exist('time_diff','var')==0
    time_diff=15;
end

if exist('qf_lower','var')==0
    qf_lower=3;
end

if exist('qf_upper','var')==0
    qf_upper=5;
end

if exist('rho','var')==0
    rho=0.028;
end

if exist('Rg','var')==0
    Rg=0.99;
end

if exist('min_null_lam','var')==0
    min_null_lam=900;
end

if exist('max_null_lam','var')==0
    max_null_lam=950;
end

%% parse out relevant header information
for ii=1:length(files)
    [data,delimiterOut]=importdata(fullfile(data_folder,files{ii}),' ',25);  %25 corresponds to number of header lines
    for jj=1:length(data.textdata)
        fn=strsplit(data.textdata{jj});
        if strcmp(fn{1},'time=')==1
            dt_ind(jj)=1;
        else
            dt_ind(jj)=0;
        end
        if strcmp(fn{1},'latitude=')==1
            lat_ind(jj)=1;
        else
            lat_ind(jj)=0;
        end
        if strcmp(fn{1},'longitude=')==1
            lon_ind(jj)=1;
        else
            lon_ind(jj)=0;
        end
    end
    ind=find(dt_ind==1);
    fn=strsplit(data.textdata{ind});
    dt(ii)=datetime([fn{5} ' ' fn{6} ' ' fn{7}],'InputFormat','MM/dd/uuuu hh:mm:ss aa'); %field 5=date, field 6=time, field 7=AM/PM
    ind=find(lat_ind==1);
    fn=strsplit(data.textdata{ind},',');
    lat_min=cell2mat(cellstr(fn{2})); %format is NNMM.DDMMA, with A indicating hemisphere; value appears to be fixed at 8 digits (4 decimal minutes)
    if length(lat_min)==0
        lat(ii)=NaN;
        lon(ii)=NaN;
        continue;
    end
    junk=str2double(lat_min(1:3))+str2double(lat_min(4:10))/60; %convert from minutes to decimal coordinates
    if strcmp(lat_min(end),'S')==1
        lat(ii)=-junk;
    else
        lat(ii)=junk;
    end
    ind=find(lon_ind==1);
    fn=strsplit(data.textdata{ind},',');
    lon_min=cell2mat(cellstr(fn{2}));  %format is NNNMM.DDMMA, with A indicating hemisphere; value appears to be fixed at 9 digits (4 decimal minutes)
    junk=str2double(lon_min(1:4))+str2double(lon_min(5:11))/60; %convert from minutes to decimal coordinates
    if strcmp(lon_min(end),'W')==1
        lon(ii)=-junk;
    else
        lon(ii)=junk;
    end
end


%identify stations

clear ff
cnt=1;
stn_cnt=0;

while cnt < length(lat)
    stn_cnt=stn_cnt+1;
    ind=find(dt<dt(cnt)+minutes(time_diff) & dt>dt(cnt)-minutes(time_diff));
    stn_lam=wavelength(:,ind);
    stn_rad=radiance(:,ind);
    for jj=1:length(ind)
        stn_files{jj}=files{ind(jj)};
    end
    svc_output(stn_cnt).latitude=nanmean(lat(ind));
    svc_output(stn_cnt).longitude=nanmean(lon(ind));
    svc_output(stn_cnt).datetime=datetime(datevec(nanmedian(datenum(dt(ind)))));
    for jj=1:length(ind)
        lind=find(wavelength(:,ind(jj)) >= 450 & wavelength(:,ind(jj)) <= 500);
        ff(jj)=nanmean(radiance(lind,ind(jj)));
    end
    [idx,schrodinger]=kmeans(ff',3);
    
    %check for quality of data by counting spectra
    ind1=find(idx==1);
    ind2=find(idx==2);
    ind3=find(idx==3);
    
    if length(ind1) < qf_lower | length(ind2) < qf_lower | length(ind3) < qf_lower | length(ind1) > qf_upper | length(ind2) > qf_upper | length(ind3) > qf_upper
        svc_output(stn_cnt).quality_flag=1;
    else
        svc_output(stn_cnt).quality_flag=0;
    end
    
    %find water spectra
    nind=find(schrodinger==min(schrodinger));
    wtr=find(idx==nind);

    lam_ind=find(wavelength(:,wtr(1)) >= 400 & wavelength(:,wtr(1)) <=600);
    
    tst=nanmedian(stn_rad(lam_ind,wtr),2);
    dev=nanstd(stn_rad(lam_ind,wtr),1,2);
    lwr=tst-dev;
    upr=tst+dev;
    
    for ii=1:length(wtr)
        if median(stn_rad(lam_ind,wtr(ii))) < median(lwr)
            igood(ii)=0;
        elseif median(stn_rad(lam_ind,wtr(ii))) > median(upr)
            igood(ii)=0;
        else
            igood(ii)=1;
        end
    end
    
    igood=find(igood==1);
    wtr=wtr(igood);
 
    svc_output(stn_cnt).water_lam=stn_lam(:,wtr);
    svc_output(stn_cnt).water_rad=stn_rad(:,wtr);
    for jj=1:length(wtr)
        svc_output(stn_cnt).water_files{jj}=stn_files{wtr(jj)};
    end
    %find sky spectra
    nind=find(schrodinger==median(schrodinger));
    sky=find(idx==nind);

    lam_ind=find(wavelength(:,sky(1)) >= 400 & wavelength(:,sky(1)) <=600);

    tst=nanmedian(stn_rad(lam_ind,sky),2);
    dev=nanstd(stn_rad(lam_ind,sky),1,2);
    lwr=tst-dev;
    upr=tst+dev;
    
    for ii=1:length(sky)
        if median(stn_rad(lam_ind,sky(ii))) < median(lwr)
            igood(ii)=0;
        elseif median(stn_rad(lam_ind,sky(ii))) > median(upr)
            igood(ii)=0;
        else
            igood(ii)=1;
        end
    end
    
    igood=find(igood==1);
    sky=sky(igood);
    
    svc_output(stn_cnt).sky_lam=stn_lam(:,sky);
    svc_output(stn_cnt).sky_rad=stn_rad(:,sky);
    for jj=1:length(sky)
        svc_output(stn_cnt).sky_files{jj}=stn_files{sky(jj)};
    end
    %find reference spectra
    nind=find(schrodinger==max(schrodinger));
    ref=find(idx==nind);

    lam_ind=find(wavelength(:,ref(1)) >= 400 & wavelength(:,ref(1)) <=600);

    tst=nanmedian(stn_rad(lam_ind,ref),2);
    dev=nanstd(stn_rad(lam_ind,ref),1,2);
    lwr=tst-dev;
    upr=tst+dev;
    
    for ii=1:length(ref)
        if median(stn_rad(lam_ind,ref(ii))) < median(lwr)
            igood(ii)=0;
        elseif median(stn_rad(lam_ind,ref(ii))) > median(upr)
            igood(ii)=0;
        else
            igood(ii)=1;
        end
    end
    
    igood=find(igood==1);
    ref=ref(igood);
    
    svc_output(stn_cnt).ref_lam=stn_lam(:,ref);
    svc_output(stn_cnt).ref_rad=stn_rad(:,ref);
    for jj=1:length(ref)
        svc_output(stn_cnt).ref_files{jj}=stn_files{ref(jj)};
    end
    cnt=cnt+length(ind);
    clear ff
    
    %calculate Rrs
    
    svc_output(stn_cnt).Rrs_lam=svc_output(stn_cnt).water_lam(:,1);
    svc_output(stn_cnt).Rrs=(nanmean(svc_output(stn_cnt).water_rad,2)-(rho.*nanmean(svc_output(stn_cnt).sky_rad,2)))./((pi.*nanmean(svc_output(stn_cnt).ref_rad,2))*(1/Rg));
    ind=find(svc_output(stn_cnt).water_lam(:,1) >= min_null_lam & svc_output(stn_cnt).water_lam(:,1) <= max_null_lam);
    svc_output(stn_cnt).Rrs_offset=(svc_output(stn_cnt).Rrs-nanmean(svc_output(stn_cnt).Rrs(ind)));

end



