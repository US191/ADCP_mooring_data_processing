function [uintfilt,vintfilt,inttim]=adcp_filt_sub(data,ui,vi,intdepvec,nhours)

sf=1/((data.time(2)-data.time(1))*24);

for iidep = 1:length(intdepvec)
    
  adcpui = ui(iidep,:);
  adcpvi = vi(iidep,:);
  valid = find(~isnan(adcpui));
  invalid = find(isnan(adcpui));
 
  if length(valid)>1 && length(valid(1):valid(end))>240; % length(invalid)>900 & 

      % Horizontal interpolation to fill gaps
      timval = data.time(valid);
      adcpui = interp1(timval,adcpui(valid),data.time);
      adcpvi = interp1(timval,adcpvi(valid),data.time);
      sum(isnan(adcpui));
      sum(isnan(adcpvi));
      % filtering
      adcpui(valid(1):valid(end)) = mfilter(adcpui(valid(1):valid(end)),sf,1/nhours,0,2*nhours,1);
      adcpvi(valid(1):valid(end)) = mfilter(adcpvi(valid(1):valid(end)),sf,1/nhours,0,2*nhours,1);
      sum(isnan(adcpui));
      sum(isnan(adcpvi));
%       adcpui(invalid) = NaN;
%       adcpvi(invalid) = NaN;
      sum(isnan(adcpui));
      uifilt(iidep,1:length(adcpui)) = adcpui;
      vifilt(iidep,1:length(adcpvi)) = adcpvi;
      % subsampling
      inttim = [ceil(data.time(1)):0.5:floor(data.time(end))];
      uintfilt(iidep,1:length(inttim)) = interp1(data.time,transpose(uifilt(iidep,:)),inttim);
      vintfilt(iidep,1:length(inttim)) = interp1(data.time,transpose(vifilt(iidep,:)),inttim);
      inttim = inttim;  
   
  elseif length(valid)<=1 || length(valid(1):valid(end))<=240
      % subsampling
      uifilt(iidep,1:length(adcpui)) = nan;
      vifilt(iidep,1:length(adcpvi)) = nan;
      inttim = [ceil(data.time(1)):0.5:floor(data.time(end))];
      uintfilt(iidep,1:length(inttim)) =nan;
      vintfilt(iidep,1:length(inttim)) =nan;
      
  else
      % Horizontal interpolation to fill gaps
      timval = data.time(valid);
      adcpui = interp1(timval,adcpui(valid),data.time);
      adcpvi = interp1(timval,adcpvi(valid),data.time); 
      % filtering
      adcpui = mfilter(adcpui,sf,1/nhours,0,2*nhours,1);
      adcpvi = mfilter(adcpvi,sf,1/nhours,0,2*nhours,1);
%       adcpui(invalid) = NaN;
%       adcpvi(invalid) = NaN;
      uifilt(iidep,1:length(adcpui)) = adcpui;
      vifilt(iidep,1:length(adcpvi)) = adcpvi;
      %subsampling
      inttim = [ceil(data.time(1)):0.5:floor(data.time(end))];
      uintfilt(iidep,1:length(inttim)) = interp1(timval,transpose(uifilt(iidep,valid)),inttim);
      vintfilt(iidep,1:length(inttim)) = interp1(timval,transpose(vifilt(iidep,valid)),inttim);
      inttim = inttim;
      
  end
      
end

hf=figure('position', [0, 0, 1400, 1000]);
subplot 121;
imagesc(data.time,intdepvec,ui,[-1 1]);
%set(gca,'YLim',[min(intdepvec)-10 max(intdepvec)+10]);
gregtick;
ylabel('Bins');
title('U field - data raw');
subplot 122;
imagesc(data.time,intdepvec,uintfilt,[-1 1]);
%set(gca,'YLim',[min(intdepvec)-10 max(intdepvec)+10]);
gregtick;
ylabel('Bins');
title('U field - data interpolated, filtered and subsampled');

hf=figure('position', [0, 0, 1400, 1000]);
subplot 121;
imagesc(data.time,intdepvec,vi,[-1 1]);
%set(gca,'YLim',[min(intdepvec)-10 max(intdepvec)+10])
gregtick;
ylabel('Bins');
title('V field - data raw');
subplot 122;
imagesc(data.time,intdepvec,vintfilt,[-1 1]);
%set(gca,'YLim',[min(intdepvec)-10 max(intdepvec)+10]);
gregtick;
ylabel('Bins');
title('V field - data interpolated, filtered and subsampled');
