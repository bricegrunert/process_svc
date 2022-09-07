function [qwip, abs_qwip] = qwip_score(rrs,rrs_wave)
% ------------
% qwip_score.m
% ------------
% Written by K. J. Turner, Tzortziou Bio-Optics Lab, City College of New
% York
% Last updated: June 8, 2022
% -------------------------------------------------------------------------
% This function calculates the Quality Water Index Polynomial (QWIP) score
% to assess the quality of hyperspectral remote sensing reflectance data
% based on spectral shape. See Dierssen et al. (2022) for more details.

% *** NOTE: This code is only applicable to hyperspectral data with 1-nm
% resolution spanning at a minimum from 400-700 nm!! See code provided in
% Dierssen et al. (2022) for application to multi-spectral satellite
% sensors. ***
% -------------------------------------------------------------------------
%   INPUTS:   
%   -------
%       rrs = remote sensing reflectance data in units of sr^-1 (vector or
%             matrix with wavelengths as rows and spectra as columns)
%
%       rrs_wave = wavelengths corresponding to rrs in units of nm (vector
%                  or matrix with same dimensions as rrs)
%
%   OUTPUTS:  
%   --------
%       qwip = QWIP score each input spectra
%
%       abs_qwip = absolute value of the QWIP score for each input spectra
% -------------------------------------------------------------------------
%   REFERENCE:
%   ----------
%   Dierssen, H. M., Vandermeulen, R. A., Barnes, B. B., Castagna, A., 
%   Knaeps, E., & Vanhellemont, Q. (2022). QWIP: A Quantitative Metric for 
%   Quality Control of Aquatic Reflectance Spectral Shape Using the Apparent 
%   Visible Wavelength. Frontiers in Remote Sensing. 3:869611. 
%   DOI: 10.3389/frsen.2022.869611.
% -------------------------------------------------------------------------

if exist('rrs','var')==0
    error('No input rrs provided');
end

if exist('rrs_wave','var')==0
    error('No input rrs wavelengths provided');
end

if size(rrs,1) ~= size(rrs_wave,1) || size(rrs,2) ~= size(rrs_wave,2)
    error('Size of rrs and rrs_wave are not the same');
end

% Extract only 400-700nm
% ----------------------
wave_min = 400;
wave_max = 700;

for ii = 1:size(rrs,2)
    wave(:,ii) = rrs_wave(rrs_wave(:,ii)>=wave_min & rrs_wave(:,ii)<=wave_max,ii);
    rrs2(:,ii) = rrs(rrs_wave(:,ii)>=wave_min & rrs_wave(:,ii)<=wave_max,ii);
end

% Calculate the Apparent Visible Wavelength (AVW)
% -----------------------------------------------
for ii = 1:size(rrs2,2)
    avw(:,ii) = sum(rrs2(:,ii))./sum(rrs2(:,ii)./wave(:,ii));
end

% Calculate the Quality Control Index (QCI)
% (Rrs_665-Rrs_490)/((Rrs_665+Rrs_490)
% -----------------------------------------
for ii = 1:size(rrs2,2)
    qci(:,ii) = (rrs2(wave(:,ii)==670,ii) - rrs2(wave(:,ii)==490,ii))./(rrs2(wave(:,ii)==670,ii) + rrs2(wave(:,ii)==490,ii));;
end

% Calculate the predicted QCI based on the polynomial provided in Dierssen
% et al. (2022) and the QWIP score as the difference bewteen the predicted
% and calculate QCI values
% ------------------------------------------------------------------------
poly_avw = 400:630;
p = [-8.399885e-09,1.715532e-05,-1.301670e-02,4.357838,-5.449532e02];
fit1 = polyval(p,poly_avw);
qci_pred = (p(1)*avw.^4 + p(2)*avw.^3 + p(3)*avw.^2 + p(4)*avw.^1 + p(5));

% Calculate outputs
% -----------------
qwip = qci_pred-qci;
abs_qwip = abs(qwip);

end 
