%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% merge_PIRATA_ADCP_10W.m
% ------------------------------------
% Merge ADCP datasets from 2001 to 2017
% -------------------------------
% Author : Jérémie HABASQUE - IRD
% -------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;

% path
addpath('.\moored_adcp_proc');

fpath = 'C:\Users\jhabasqu\Desktop\PIRATA\ADCP_mooring_data\10W-0N\merge_data\';
mooring.name='10W0N';

% niveaux d'affichages des champs u et v
niv_u = (-1.5:0.1:1.5);
niv_v = (-0.5:0.1:0.5);

% subsampling on a regular 24-hour grid
% 2001,2003 and 2004 datasets are processed with a daily resolution
step_subsampling = 1; 

% Read data

% 2001 (on le lit pas car periode de 2 mois..)
load('C:\Users\jhabasqu\Desktop\PIRATA\ADCP_mooring_data\10W-0N\2001-2002\UV_day_10w_int.mat');
d0=datenum(2001,12,11,0,0,0);
d1=datenum(2002,02,23,0,0,0); 
adcp_2001_time=d0:step_subsampling:d1;
adcp_2001_time = adcp_2001_time';
adcp_2001_u = uvmoy(:,1:75)/100;
adcp_2001_v = uvmoy(:,76:150)/100;
adcp_2001_bin_length = 4;
adcp_2001_z = 0:adcp_2001_bin_length:104;%Pas certain...
adcp_2001_z = adcp_2001_z';%Pas certain...  
[YY,MM,DD,hh,mm,ss] = datevec(adcp_2001_time);
adcp_2001_time=julian(YY,MM,DD,hh,mm,ss);
 
% 2003 
load('C:\Users\jhabasqu\Desktop\PIRATA\ADCP_mooring_data\10W-0N\2003-2004\UV_moy_10w.mat');
d0=datenum(2003,05,07,0,0,0);
d1=datenum(2004,02,03,0,0,0); 
adcp_2003_time=d0:step_subsampling:d1;
adcp_2003_time = adcp_2003_time';
adcp_2003_u = uvmoy(1:39,:)/100;
adcp_2003_v = uvmoy(40:78,:)/100;
adcp_2003_bin_length = 8;
adcp_2003_z = 26:adcp_2003_bin_length:332;%Pas certain...
adcp_2003_z = adcp_2003_z';%Pas certain...  
[YY,MM,DD,hh,mm,ss] = datevec(adcp_2003_time);
adcp_2003_time=julian(YY,MM,DD,hh,mm,ss);

% 2004
load('C:\Users\jhabasqu\Desktop\PIRATA\ADCP_mooring_data\10W-0N\2004-2005\UV_day_10w_int10.mat');
d0=datenum(2004,02,06,0,0,0);
d1=datenum(2005,06,17,0,0,0); 
adcp_2004_time=d0:step_subsampling:d1;
adcp_2004_time = adcp_2004_time';
adcp_2004_u=uvmoy(:,1:498)/100;
adcp_2004_v=uvmoy(:,499:996)/100;
adcp_2004_bin_length = 10;
adcp_2004_z = 0:adcp_2004_bin_length:301;%Pas certain...
adcp_2004_z = adcp_2004_z';%Pas certain...
[YY,MM,DD,hh,mm,ss] = datevec(adcp_2004_time);
adcp_2004_time=julian(YY,MM,DD,hh,mm,ss);

%2006-2008
load('C:\Users\jhabasqu\Desktop\PIRATA\ADCP_mooring_data\10W-0N\2006-2008\PIRATA_FR18_10W_daily.mat');
d0=datenum(2006,06,26,00,00,0);
d1=datenum(2008,09,27,00,00,0); 
adcp_2006_time=d0:step_subsampling:d1;
adcp_2006_time = adcp_2006_time';
[YY,MM,DD,hh,mm,ss] = datevec(adcp_2006_time);
adcp_2006_time=julian(YY,MM,DD,hh,mm,ss);
adcp_2006_u = u'/100;
adcp_2006_v = v'/100;
adcp_2006_z = Z'; 

% 2011
load('C:\Users\jhabasqu\Desktop\PIRATA\ADCP_mooring_data\10W-0N\2010-2011\FR22-10W0N_UP_DOWN_int_filt_sub.mat');
adcp_2011_time = data.inttim;
adcp_2011_u = data.u_final; 
adcp_2011_v = data.v_final; 
adcp_2011_z = data.z_final';
adcp_2011_bin_length = raw.config.cell;

% 2014
load('C:\Users\jhabasqu\Desktop\PIRATA\ADCP_mooring_data\10W-0N\2012-2013\FR24-10W_15258_instr_01_int_filt_sub.mat');
adcp_2014_time = data.inttim;
adcp_2014_u = data.uintfilt;
adcp_2014_v = data.vintfilt;
adcp_2014_z = data.Z';
adcp_2014_bin_length = raw.config.cell;

% 2015
load('C:\Users\jhabasqu\Desktop\PIRATA\ADCP_mooring_data\10W-0N\2014-2015\FR25-10W_15258_instr_01_int_filt_sub.mat');
adcp_2015_time = data.inttim;
adcp_2015_u = data.uintfilt;
adcp_2015_v = data.vintfilt;
adcp_2015_z = data.Z';
adcp_2015_bin_length = raw.config.cell;

% 2017
load('C:\Users\jhabasqu\Desktop\PIRATA\ADCP_mooring_data\10W-0N\2015-2017\10W0N_15258_instr_01_int_filt_sub.mat');
adcp_2017_time = data.inttim;
adcp_2017_u = data.uintfilt;
adcp_2017_v = data.vintfilt;
adcp_2017_z = data.Z';
adcp_2017_bin_length = raw.config.cell;

%% Interpolate data on a regular vertical grid

%Z = fliplr(min_depth:max_bin_length:max_depth);
Z = fliplr(0:5:350); % on se base sur la maille des donnees TACE etendue a 350m
Zmax = max(Z);

%interpolation for each timestep for 2001 data
u_interp_2001 = NaN(length(Z),length(adcp_2001_time));
v_interp_2001 = NaN(length(Z),length(adcp_2001_time));
for i=1:length(adcp_2001_time)    
     ind_ok = find(~isnan(adcp_2001_u(:,i)));
     u_interp_2001(:,i) = interp1(adcp_2001_z(ind_ok),adcp_2001_u(ind_ok,i),Z);
     v_interp_2001(:,i) = interp1(adcp_2001_z(ind_ok),adcp_2001_v(ind_ok,i),Z);    
end
%interpolation for each timestep for 2003 data
u_interp_2003 = NaN(length(Z),length(adcp_2003_time));
v_interp_2003 = NaN(length(Z),length(adcp_2003_time));
for i=1:length(adcp_2003_time)    
     ind_ok = find(~isnan(adcp_2003_u(:,i)));
     u_interp_2003(:,i) = interp1(adcp_2003_z(ind_ok),adcp_2003_u(ind_ok,i),Z);
     v_interp_2003(:,i) = interp1(adcp_2003_z(ind_ok),adcp_2003_v(ind_ok,i),Z);    
end

%interpolation for each timestep for 2004 data
u_interp_2004 = NaN(length(Z),length(adcp_2004_time));
v_interp_2004 = NaN(length(Z),length(adcp_2004_time));
for i=1:length(adcp_2004_time)    
     ind_ok = find(~isnan(adcp_2004_u(:,i)));
     u_interp_2004(:,i) = interp1(adcp_2004_z(ind_ok),adcp_2004_u(ind_ok,i),Z);
     v_interp_2004(:,i) = interp1(adcp_2004_z(ind_ok),adcp_2004_v(ind_ok,i),Z);    
end

%interpolation for each timestep for 2006 data
u_interp_2006 = NaN(length(Z),length(adcp_2006_time));
v_interp_2006 = NaN(length(Z),length(adcp_2006_time));
for i=1:length(adcp_2006_time)    
     ind_ok = find(~isnan(adcp_2006_u(:,i)));
     u_interp_2006(:,i) = interp1(adcp_2006_z(ind_ok),adcp_2006_u(ind_ok,i),Z);
     v_interp_2006(:,i) = interp1(adcp_2006_z(ind_ok),adcp_2006_v(ind_ok,i),Z);    
end

%interpolation for each timestep for 2011 data
u_interp_2011 = NaN(length(Z),length(adcp_2011_time));
v_interp_2011 = NaN(length(Z),length(adcp_2011_time));
for i=1:length(adcp_2011_time)    
     ind_ok = find(~isnan(adcp_2011_u(:,i)));
     u_interp_2011(:,i) = interp1(adcp_2011_z(ind_ok),adcp_2011_u(ind_ok,i),Z);
     v_interp_2011(:,i) = interp1(adcp_2011_z(ind_ok),adcp_2011_v(ind_ok,i),Z);  
end

%interpolation for each timestep for 2014 data
u_interp_2014 = NaN(length(Z),length(adcp_2014_time));
v_interp_2014 = NaN(length(Z),length(adcp_2014_time));
for i=1:length(adcp_2014_time)
    ind_ok = find(~isnan(adcp_2014_u(:,i)));
    u_interp_2014(:,i) = interp1(adcp_2014_z(ind_ok),adcp_2014_u(ind_ok,i),Z);
    v_interp_2014(:,i) = interp1(adcp_2014_z(ind_ok),adcp_2014_v(ind_ok,i),Z);
end

%interpolation for each timestep for 2015 data
u_interp_2015 = NaN(length(Z),length(adcp_2015_time));
v_interp_2015 = NaN(length(Z),length(adcp_2015_time));
for i=1:length(adcp_2015_time)
 ind_ok = find(~isnan(adcp_2015_u(:,i)));
    u_interp_2015(:,i) = interp1(adcp_2015_z(ind_ok),adcp_2015_u(ind_ok,i),Z);
    v_interp_2015(:,i) = interp1(adcp_2015_z(ind_ok),adcp_2015_v(ind_ok,i),Z);
end
 
%interpolation for each timestep for 2017 data
u_interp_2017 = NaN(length(Z),length(adcp_2017_time));
v_interp_2017 = NaN(length(Z),length(adcp_2017_time));
for i=1:length(adcp_2017_time)
 ind_ok = find(~isnan(adcp_2017_u(:,i)));
    u_interp_2017(:,i) = interp1(adcp_2017_z(ind_ok),adcp_2017_u(ind_ok,i),Z);
    v_interp_2017(:,i) = interp1(adcp_2017_z(ind_ok),adcp_2017_v(ind_ok,i),Z);
end

%% combine all data
all_time = [adcp_2001_time' adcp_2003_time' adcp_2004_time' adcp_2006_time' adcp_2011_time adcp_2014_time adcp_2015_time adcp_2017_time];
%check sampling interval
% figure;
% hist(diff(all_time),100); xlim([0, 2]);

all_u_interp = [u_interp_2001 u_interp_2003 u_interp_2004 u_interp_2006 u_interp_2011 u_interp_2014 u_interp_2015 u_interp_2017];
all_v_interp = [v_interp_2001 v_interp_2003 v_interp_2004 v_interp_2006 v_interp_2011 v_interp_2014 v_interp_2015 v_interp_2017];

% create a continuous series of daily data, ranging from min(d) to max(d)
ADCP_final.time = ceil(min(all_time)):step_subsampling:floor(max(all_time));
ADCP_final.depth = Z;
ADCP_final.u = NaN(length(ADCP_final.depth),length(ADCP_final.time));
ADCP_final.v = NaN(length(ADCP_final.depth),length(ADCP_final.time));

for i_time = 1:length(ADCP_final.time)
    for j_time = 1:length(all_time)
        if ADCP_final.time(i_time) == all_time(j_time)
            ADCP_final.u(:,i_time)=all_u_interp(:,j_time);
            ADCP_final.v(:,i_time)=all_v_interp(:,j_time);
        end
    end
end

% save global data
save([fpath, 'ADCP_10W0N_2001_2017_int_filt_sub.mat'],'ADCP_final');

%% FIGURES

hf=figure('position', [0, 0, 1400, 1000]);
%u
subplot(2,1,1);
[C,h] = contourf(ADCP_final.time,Z,ADCP_final.u,niv_u);
set(h,'LineColor','none');
caxis(niv_u([1 end]));
h=colorbar;
ylabel(h,'U [m s^-^1]');
set(gca,'ydir', 'reverse');
ylim([0 350]);
ylabel('Depth (m)');
%change figure label in HH:MM
gregtick
title({['10°W 0°N - ZONAL VELOCITY']});

%v
subplot(2,1,2);
[C,h] = contourf(ADCP_final.time,Z,ADCP_final.v,niv_v);
set(h,'LineColor','none');
caxis(niv_v([1 end]));
h=colorbar;
ylabel(h,'V [m s^-^1]');
set(gca,'ydir', 'reverse');
ylim([0 350]);
ylabel('Depth (m)');
%change figure label in HH:MM
gregtick
title({['10°W 0°N - MERIDIONAL VELOCITY']});

graph_name = [fpath, 'ADCP_10W0N_2001_2017_U_V_daily'];
set(hf,'Units','Inches');
pos = get(hf,'Position');
set(hf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hf,graph_name,'-dpdf','-r300');

% Histogramme des valeurs U et V

hf=figure('position', [0, 0, 1400, 1000]);
subplot(1,2,1); hist(ADCP_final.u(:),100); xlabel('U [m s^-^1]');
subplot(1,2,2); hist(ADCP_final.v(:),100); xlabel('V [m s^-^1]');

graph_name = [fpath, 'ADCP_U_V_10W_daily_histo'];
set(hf,'Units','Inches');
pos = get(hf,'Position');
set(hf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hf,graph_name,'-dpdf','-r300');


%%  Write netcdf file                           
ncid=netcdf.create([fpath,'ADCP_10W0N_2001_2017_1d.nc'],'NC_WRITE');
 
%create dimension
dimidt = netcdf.defDim(ncid,'time',length(ADCP_final.time));
dimidz = netcdf.defDim(ncid,'depth',length(ADCP_final.depth));
%Define IDs for the dimension variables (pressure,time,latitude,...)
time_ID=netcdf.defVar(ncid,'time','double',[dimidt]);
depth_ID=netcdf.defVar(ncid,'depth','double',[dimidz]);
%Define the main variable ()
u_ID = netcdf.defVar(ncid,'u','double',[dimidt dimidz]);
v_ID = netcdf.defVar(ncid,'v','double',[dimidt dimidz]);
%We are done defining the NetCdf
netcdf.endDef(ncid);
%Then store the dimension variables in
netcdf.putVar(ncid,time_ID,ADCP_final.time);
netcdf.putVar(ncid,depth_ID,ADCP_final.depth);  
%Then store my main variable
netcdf.putVar(ncid,u_ID,ADCP_final.u);
netcdf.putVar(ncid,v_ID,ADCP_final.v);
%We're done, close the netcdf
netcdf.close(ncid);
