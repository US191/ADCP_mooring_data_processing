function [uintfilt,vintfilt,uifilt,vifilt,inttim,utid_baro,vtid_baro]=adcp_filt_sub(data,ui,vi,intdepvec,nhours)

sf = 1/((data.time(2)-data.time(1))*24);

for iidep = 1:length(intdepvec)
    
  adcpui  = ui(iidep,:);
  adcpvi  = vi(iidep,:);
  valid   = find(~isnan(adcpui));
  invalid = find(isnan(adcpui));
 
  if length(valid)>1 && length(valid(1):valid(end))>240 % length(invalid)>900 & 

      % Horizontal interpolation to fill gaps
      timval                           = data.time(valid);
      adcpui                           = interp1(timval,adcpui(valid),data.time);
      adcpvi                           = interp1(timval,adcpvi(valid),data.time);
      sum(isnan(adcpui));
      sum(isnan(adcpvi));
      % filtering
      adcpui(valid(1):valid(end))      = mfilter(adcpui(valid(1):valid(end)),sf,1/nhours,0,2*nhours,1);
      adcpvi(valid(1):valid(end))      = mfilter(adcpvi(valid(1):valid(end)),sf,1/nhours,0,2*nhours,1);
      %[adcpui(valid(1):valid(end)) ] = fft_filter(adcpui(valid(1):valid(end)),2,[24 Inf]);
      %[adcpvi(valid(1):valid(end))]  = fft_filter(adcpvi(valid(1):valid(end)),2,[24 Inf]);
      sum(isnan(adcpui));
      sum(isnan(adcpvi));
%       adcpui(invalid) = NaN;
%       adcpvi(invalid) = NaN;
      sum(isnan(adcpui));
      uifilt(iidep,1:length(adcpui))   = adcpui;
      vifilt(iidep,1:length(adcpvi))   = adcpvi;
      % subsampling
      inttim                           = [ceil(data.time(1)):0.25:floor(data.time(end))];
      uintfilt(iidep,1:length(inttim)) = interp1(data.time,transpose(uifilt(iidep,:)),inttim);
      vintfilt(iidep,1:length(inttim)) = interp1(data.time,transpose(vifilt(iidep,:)),inttim);
      inttim = inttim;  
   
  elseif length(valid)<=1 || length(valid(1):valid(end))<=240
      % subsampling
      uifilt(iidep,1:length(adcpui))   = NaN;
      vifilt(iidep,1:length(adcpvi))   = NaN;
      inttim                           = [ceil(data.time(1)):0.25:floor(data.time(end))];
      uintfilt(iidep,1:length(inttim)) = NaN;
      vintfilt(iidep,1:length(inttim)) = NaN;
      
  else
      % Horizontal interpolation to fill gaps
      timval                           = data.time(valid);
      adcpui                           = interp1(timval,adcpui(valid),data.time);
      adcpvi                           = interp1(timval,adcpvi(valid),data.time); 
      %adcpui                           = mfilter(adcpui,sf,1/nhours,0,2*nhours,1);
      %adcpvi                           = mfilter(adcpvi,sf,1/nhours,0,2*nhours,1);
%       adcpui(invalid) = NaN;
%       adcpvi(invalid) = NaN;
      uifilt(iidep,1:length(adcpui))   = adcpui;
      vifilt(iidep,1:length(adcpvi))   = adcpvi;
      %subsampling
      inttim                           = [ceil(data.time(1)):0.25:floor(data.time(end))];
      uintfilt(iidep,1:length(inttim)) = interp1(timval,transpose(uifilt(iidep,valid)),inttim);
      vintfilt(iidep,1:length(inttim)) = interp1(timval,transpose(vifilt(iidep,valid)),inttim);
      inttim                           = inttim;
      
  end
      
end

% [uintfilt] = fft_filter(uintfilt,2,[40 Inf]);
% [vintfilt] = fft_filter(vintfilt,2,[40 Inf]);
% uintfilt   = real(uintfilt);
% vintfilt   = real(vintfilt);

%% Calc tide
addpath(genpath('C:\Users\proussel\Documents\outils\TOOLS\TMD'))
addpath(genpath('C:\Users\proussel\Documents\outils\TOOLS\TIDES\tpxo9_atlas'))
addpath(genpath('C:\Users\proussel\Documents\outils\TOOLS\TIDES\tpxo8_atlas'))
% u component
utid      = tmd_tide_pred('C:\Users\proussel\Documents\outils\TOOLS\TIDES\tpxo8_atlas_compact\Model_tpxo8.v1',datenum(julian(inttim')),data.lat,data.lon,'u')/100;
%utid      = tmd_tide_pred('C:\Users\proussel\Documents\outils\TOOLS\TIDES\tpxo9.1\Model_tpxo9.v1',datenum(julian(inttim')),data.lat,data.lon,'u')/100;
%utid      = tmd_tide_pred('C:\Users\proussel\Documents\outils\TOOLS\TMD\DATA\Model_tpxo7.2',datenum(julian(inttim')),data.lat,data.lon,'u')/100;
utid_baro = utid'.*ones(length(uintfilt(:,1)),1);
% v component
vtid      = tmd_tide_pred('C:\Users\proussel\Documents\outils\TOOLS\TIDES\tpxo8_atlas_compact\Model_tpxo8.v1',datenum(julian(inttim')),data.lat,data.lon,'v')/100;
%vtid      = tmd_tide_pred('C:\Users\proussel\Documents\outils\TOOLS\TIDES\tpxo9.1\Model_tpxo9.v1',datenum(julian(inttim')),data.lat,data.lon,'v')/100;
%vtid      = tmd_tide_pred('C:\Users\proussel\Documents\outils\TOOLS\TMD\DATA\Model_tpxo7.2',datenum(julian(inttim')),data.lat,data.lon,'v')/100;
vtid_baro = vtid'.*ones(length(vintfilt(:,1)),1);
utid_baro=0;
vtid_baro=0;

%% Correct tide
% uintfilt = uintfilt-utid_baro;
% vintfilt = vintfilt-vtid_baro;

hf=figure('position', [0, 0, 1400, 1000]);
subplot 121;
imagesc(data.time,intdepvec,ui(intdepvec,:),[-1 1]);
%set(gca,'YLim',[min(intdepvec)-10 max(intdepvec)+10]);
gregtick;
ylabel('Bins');
title('U field - data raw');
subplot 122;
imagesc(inttim,intdepvec,uintfilt,[-1 1]);
%set(gca,'YLim',[min(intdepvec)-10 max(intdepvec)+10]);
gregtick;
ylabel('Bins');
title('U field - data interpolated, filtered and subsampled');

hf=figure('position', [0, 0, 1400, 1000]);
subplot 121;
imagesc(data.time,intdepvec,vi(intdepvec,:),[-1 1]);
%set(gca,'YLim',[min(intdepvec)-10 max(intdepvec)+10])
gregtick;
ylabel('Bins');
title('V field - data raw');
subplot 122;
imagesc(inttim,intdepvec,vintfilt,[-1 1]);
%set(gca,'YLim',[min(intdepvec)-10 max(intdepvec)+10]);
gregtick;
ylabel('Bins');
title('V field - data interpolated, filtered and subsampled');
