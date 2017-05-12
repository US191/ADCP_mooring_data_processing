%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% template_get_adcp_data.m
% -------------------------------
% Author : J�r�mie HABASQUE - IRD
% -------------------------------
% INPUTS:
% - binary raw file with .000 extension
% OUTPUTS:
% - U and V fields interpolated on a regulard grid, filtered and subsampled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First part --------------------------------------------------------------------------------------------------------------------
close all
clear all
 
%% META information:

% Path
addpath('.\moored_adcp_proc'); % ou par exemple ('C:\Users\IRD_US_IMAGO\TRAITEMENTS\ADCP_MOUILLAGE\01_DATA_PROCESSING\moored_adcp_proc');
addpath('.\backscatter'); % (Optionnel) / ou par exemple ('C:\Users\IRD_US_IMAGO\TRAITEMENTS\ADCP_MOUILLAGE\01_DATA_PROCESSING\backscatter');

% Location rawfile
fpath = '';
rawfile='.\data_example\FR24_000.000'; % binary file with .000 extension
 
% Directory for outputs
fpath_output = '.\data_example\';

% Cruise/mooring info
cruise.name = 'Cruise Name';
mooring.name='Lat Lon'; % 0N10W par exemple
mooring.lat=00+00/60; %latitude en degr�s d�cimaux
mooring.lon=-10+00/60; %longitude en degr�s d�cimaux

% ADCP info
adcp.sn=15258;
adcp.type='150 khz Quartermaster'; % Type : �Quartermaster�, �longranger�
adcp.direction='up';        % upward-looking 'up', downward-looking 'dn'
adcp.instr_depth=280;       % nominal instrument depth
instr=1;                    % this is just for name convention and sorting of all mooring instruments

% If ADCP was not set up to correct for magnetic deviation internally
% ("EA0" code in configuration file), use http://www.ngdc.noaa.gov/geomag-web/#declination
% Magnetic deviation: Mean of deviations at time of deployment and time of recovery 

% Magnetic deviation values
magnetic_deviation_ini = 9.29;
magnetic_deviation_end = 9.05;
rot=-(magnetic_deviation_ini+magnetic_deviation_end)/2;  

% Read rawfile
fprintf('Read %s\n', rawfile);
raw=read_os3(rawfile,'all');
figure;plot(raw.pressure);set(gca,'ydir','reverse');
title('pressure sensor');ylabel('Depth(m)');xlabel('Time index');grid on; 
saveas(gcf,[fpath_output,mooring.name,'_',num2str(adcp.sn),'_instr_',num2str(instr),'_','Pressure_sensor'],'fig')

% Second part --------------------------------------------------------------------------------------------------------------------

% Determine first and last indiced when instrument was at depth (you can do this by plotting 'raw.pressure' for example           
first = 12; 
last = 17620; 

% amplitude of the bins / Correction ADCP's depth
ea = squeeze(mean(raw.amp(:,:,first:last),2));   
figure; imagesc(ea);title('Amplitude of the bins'); colorbar;
ylabel('Bins');xlabel('Time index'); 
saveas(gcf,[fpath_output,mooring.name,'_',num2str(adcp.sn),'_instr_',num2str(instr),'_','Amplitude_bins'],'fig')

% Third part --------------------------------------------------------------------------------------------------------------------

% If upward looking: range of surface bins used for instrument depth correction below!
sbins= 31:38; % here a range of bins is given which cover the surface reflection

% Exclude data with percent good below prct_good
prct_good = 20;

%% Read data
freq = raw.config.sysconfig.frequency;

u2 = squeeze(raw.vel(:,1,first:last));
v2 = squeeze(raw.vel(:,2,first:last));
w = squeeze(raw.vel(:,3,first:last));
vel_err = squeeze(raw.vel(:,4,first:last));  % the difference in vertical velocity between the two pairs of transducers 
time = raw.juliandate(first:last);
ang = [raw.pitch(first:last) raw.roll(first:last) raw.heading(first:last)]; 
soundspeed = raw.soundspeed(first:last);
T = raw.temperature(first:last);
press = raw.pressure(first:last);

nbin = raw.config.ncells;  % number of bins
bin  = 1:nbin;
blen = raw.config.cell;    % bin length
blnk = raw.config.blank;   % blank distance after transmit

dt=(time(2)-time(1))*24;   % Sampling interval in hours

[u,v]=uvrot(u2,v2,-rot);   % Correction of magnetic deviation

pg = squeeze(raw.pg(:,4,first:last));  % percent good

igap=find(pg<prct_good);           % Exclude data with percent good below prct_good
u(igap)=nan;
v(igap)=nan;
w(igap)=nan;
vel_err(igap)=nan;

%% Calculate depth of each bin:
dpt = sw_dpth(press,mooring.lat)';  % convert pressure to depth, press needs to have dimension (n x 1)
dpt1 = repmat(dpt,nbin,1);
binmat = repmat((1:nbin)',1,length(dpt1));

% If ADCP is upward-looking a depth correction can be inferred from the
% surface reflection, which is done in adcp_surface_fit
if strcmp(adcp.direction,'up')  
    [z,dpt1,offset,xnull]=adcp_surface_fit(dpt,ea,sbins,blen,blnk,nbin);
elseif strcmp(adcp.direction,'dn')
    z = dpt1+(binmat-0.5)*blen+blnk;
else
    error('Bin depth calculation: unknown direction!');
end

saveas(figure(1),[fpath_output,mooring.name,'_',num2str(adcp.sn),'_instr_',num2str(instr),'_','Hist_diff_orig-depth_recon-depth'],'fig')
saveas(figure(2),[fpath_output,mooring.name,'_',num2str(adcp.sn),'_instr_',num2str(instr),'_','Offset_depth'],'fig')
saveas(figure(3),[fpath_output,mooring.name,'_',num2str(adcp.sn),'_instr_',num2str(instr),'_','Amplitude_bins_2'],'fig')

%% Remove bad data if ADCP is looking upward
u1=u; v1=v; w1=w; vel_err1=vel_err; ea1=ea;

if strcmp(adcp.direction,'up')
    for i=1:length(time)
        sz_dpt(i)=adcp_shadowzone(dpt(i),raw.config.sysconfig.angle); % depending on the instrument depth and the beam angle the shadow zone, i.e. the depth below the surface which is contaminated by the surface reflection is determined

        iz(i)=find(z(:,i)>sz_dpt(i),1,'last');
        sbin(i)=bin(iz(i));
        
        % here a manual criterion should be hard-coded if
        % adcp_check_surface (below) shows bad velocities close to the
        % surface

        fz(i)=z(iz(i),i);    
    end

    for i=1:length(time) 
        u1(sbin(i)+1:end,i)=nan;
        v1(sbin(i)+1:end,i)=nan;
        w1(sbin(i)+1:end,i)=nan;
        vel_err1(sbin(i)+1:end,i)=nan;
        ea1(sbin(i)+1:end,i)=nan;
    end

    if 1
        bins=nmedian(sbin)-4:nmedian(sbin)+2;
        adcp_check_surface(bins,u,u1,v,v1,time,bin,z); 
        % here the closest bins below the surface are plotted that are supposed to have good velocities, if there are still bad velocities a manual criterion needs to be found
    end
end
saveas(figure(4),[fpath_output,mooring.name,'_',num2str(adcp.sn),'_instr_',num2str(instr),'_','Meridional_zonal_velocity'],'fig')

%% SAVE DATA
% More meta information
%adcp.comment='';
adcp.config=raw.config;
adcp.z_offset=offset;
adcp.ang=ang;
adcp.mag_dev=rot;

% Data structure
data.u=u1;
data.v=v1;
data.w=w1;
data.e=vel_err1;
data.ea=ea1;
data.pg=pg;
data.time=time;
data.z_bins=z;
data.depth=dpt;
data.temp=T;
data.sspd = soundspeed;

% Save
save([fpath_output, mooring.name '_' num2str(adcp.sn) '_instr_' sprintf('%02d',instr) '.mat'],'adcp','mooring','data','raw');

%% Interpolate data on a regular vertical grid
Z = fliplr(blen/2:blen:max(z(:))+blen);
Zmax = max(Z);

u_interp = NaN(length(time),length(Z));
v_interp = NaN(length(time),length(Z));
for i=1:length(time)
    % indice correspondant sur la grille finale Z
    ind = round((Zmax-z(1,i))/blen)+1;
    % filling the grid
    npts = min([length(Z)-ind+1 length(bin)]);
    u_interp(i,ind:ind+npts-1) = u1(1:npts,i);
    v_interp(i,ind:ind+npts-1) = v1(1:npts,i);
end
  
%% Horizontal interpolation, filtering and subsampling
[uintfilt,vintfilt,inttim] = adcp_filt_sub(data,u_interp',v_interp',1:length(Z),40);
saveas(figure(5),[fpath_output,mooring.name,'_',num2str(adcp.sn),'_instr_',num2str(instr),'_','data_raw_filt_subsampled_1'],'fig')
saveas(figure(6),[fpath_output,mooring.name,'_',num2str(adcp.sn),'_instr_',num2str(instr),'_','data_raw_filt_subsampled_2'],'fig')

% Save interpolated data
bin_start = 1; % bin indice where good interpolated data for the whole dataset start
data.uintfilt=uintfilt(bin_start:length(Z),:);
data.vintfilt=vintfilt(bin_start:length(Z),:);
data.Z = Z(bin_start:length(Z));
data.inttim = inttim;
save([fpath_output, mooring.name '_' num2str(adcp.sn) '_instr_' sprintf('%02d',instr) '_int_filt_sub.mat'],'adcp','mooring','data','raw');

%% Figure
niv_u = (-1.5:0.1:1.5);
niv_v = (-0.5:0.1:0.5);

hf=figure('position', [0, 0, 1400, 1000]);
%u
subplot(2,1,1);
[C,h] = contourf(inttim,Z(bin_start:length(Z)),uintfilt(bin_start:length(Z),:),niv_u); 
set(h,'LineColor','none');
caxis(niv_u([1 end]));
h=colorbar;
ylabel(h,'U [m s^-^1]');
set(gca,'ydir', 'reverse');
ylabel('Depth (m)');
ylim([0,adcp.instr_depth]);
%change figure label in HH:MM
gregtick;
title({[mooring.name, ' - MERIDIONAL VELOCITY - RDI ',num2str(freq),' kHz']});

%v
subplot(2,1,2);
[C,h] = contourf(inttim,Z(bin_start:length(Z)),vintfilt(bin_start:length(Z),:),niv_v); 
set(h,'LineColor','none');
caxis(niv_v([1 end]));
h=colorbar;
ylabel(h,'V [m s^-^1]');
set(gca,'ydir', 'reverse');
ylabel('Depth (m)');
ylim([0,adcp.instr_depth]);
%change figure label in HH:MM
gregtick;
title({[mooring.name, ' - ZONAL VELOCITY - RDI ',num2str(freq),' kHz']});

graph_name = [fpath_output, mooring.name '_U_V_int_filt_sub'];
set(hf,'Units','Inches');
pos = get(hf,'Position');
set(hf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(hf,graph_name,'-dpdf','-r300');
 
%%  Write netcdf file           
[yr_start , ~, ~] = gregorian(inttim(1));
[yr_end,  ~, ~] = gregorian(inttim(length(inttim)));

ncid=netcdf.create([fpath_output,'ADCP_',mooring.name,'_',num2str(yr_start),'_',num2str(yr_end),'_1d.nc'],'NC_WRITE');
 
%create dimension
dimidt = netcdf.defDim(ncid,'time',length(inttim));
dimidz = netcdf.defDim(ncid,'depth',length(Z));
%Define IDs for the dimension variables (pressure,time,latitude,...)
time_ID=netcdf.defVar(ncid,'time','double',dimidt);
depth_ID=netcdf.defVar(ncid,'depth','double',dimidz);
%Define the main variable ()
u_ID = netcdf.defVar(ncid,'u','double',[dimidt dimidz]);
v_ID = netcdf.defVar(ncid,'v','double',[dimidt dimidz]);
%We are done defining the NetCdf
netcdf.endDef(ncid);
%Then store the dimension variables in
netcdf.putVar(ncid,time_ID,inttim);
netcdf.putVar(ncid,depth_ID,Z);  
%Then store my main variable
netcdf.putVar(ncid,u_ID,uintfilt);
netcdf.putVar(ncid,v_ID,vintfilt);
%We're done, close the netcdf
netcdf.close(ncid);

% rmpath
% rmpath('..\moored_adcp_proc');
clear all; close all;

% -------------------------------------------------------------------------------------------