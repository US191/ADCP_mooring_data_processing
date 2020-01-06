%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% merge_PIRATA_ADCP_10W.m
% ------------------------------------
% Merge ADCP datasets from 2001 to 2019
% -------------------------------
% Author : Jérémie HABASQUE - IRD
%          Pierre Rousselot - IRD
% -------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;

% path
addpath('.\moored_adcp_proc');

fpath        = 'C:\Users\proussel\Documents\outils\ADCP\ADCP_mooring_data_processing\FR29\merge_data\';
mooring.name = '0N10W';

% niveaux d'affichages des champs u et v
niv_u = (-1.5:0.1:1.5);
niv_v = (-0.5:0.1:0.5);

% subsampling on a regular 24-hour grid
% 2001,2003 and 2004 datasets are processed with a daily resolution
step_subsampling = 1; 

% Read data

% 2001-2017
load('C:\Users\proussel\Documents\outils\ADCP\ADCP_mooring_data_processing\FR29\ADCP_10W0N_2001_2017_int_filt_sub.mat');
u     = ADCP_final.u;
v     = ADCP_final.v;
depth = ADCP_final.depth;
time  = ADCP_final.time;
clear ADCP_final

% 2019
load('C:\Users\proussel\Documents\outils\ADCP\ADCP_mooring_data_processing\FR29\0N10W_15258_instr_01_int_filt_sub.mat');
adcp_2019_time       = data.inttim;
adcp_2019_u          = data.uintfilt;
adcp_2019_v          = data.vintfilt;
adcp_2019_z          = data.Z';

%% Interpolate data on a regular vertical grid

%Z = fliplr(min_depth:max_bin_length:max_depth);
Z    = fliplr(0:5:350); % on se base sur la maille des donnees TACE etendue a 350m
Zmax = max(Z);

%interpolation for each timestep for 2019 data
u_interp_2019 = NaN(length(Z),length(adcp_2019_time));
v_interp_2019 = NaN(length(Z),length(adcp_2019_time));
for i=1:length(adcp_2019_time)
    ind_ok             = find(~isnan(adcp_2019_u(:,i)));
    u_interp_2019(:,i) = interp1(adcp_2019_z(ind_ok),adcp_2019_u(ind_ok,i),Z);
    v_interp_2019(:,i) = interp1(adcp_2019_z(ind_ok),adcp_2019_v(ind_ok,i),Z);
end

%% combine all data
all_time = [time adcp_2019_time];
%check sampling interval
% figure;
% hist(diff(all_time),100); xlim([0, 2]);

all_u_interp = [u u_interp_2019];
all_v_interp = [v v_interp_2019];

% create a continuous series of daily data, ranging from min(d) to max(d)
ADCP_final.time  = ceil(min(all_time)):step_subsampling:floor(max(all_time));
ADCP_final.depth = Z;
ADCP_final.u     = NaN(length(ADCP_final.depth),length(ADCP_final.time));
ADCP_final.v     = NaN(length(ADCP_final.depth),length(ADCP_final.time));

for i_time = 1:length(ADCP_final.time)
    for j_time = 1:length(all_time)
        if ADCP_final.time(i_time) == all_time(j_time)
            ADCP_final.u(:,i_time) = all_u_interp(:,j_time);
            ADCP_final.v(:,i_time) = all_v_interp(:,j_time);
        end
    end
end

% save global data
save([fpath, 'ADCP_10W0N_2001_2019_int_filt_sub.mat'],'ADCP_final');

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

graph_name = [fpath, 'ADCP_10W0N_2001_2019_U_V_daily'];
set(hf,'Units','Inches');
pos = get(hf,'Position');
set(hf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hf,graph_name,'-dpdf','-r300');

% % FIGURE 2
% cUaxis = [-0.8 0.8];   % For U colorbar
% cVaxis = [-0.4 0.4];   % For V colorbar
% parU = jet(64);
% parV = jet(32); 
% figure
% subplot(2,1,1)
% pcolor(ADCP_final.time,Z,ADCP_final.u); shading interp
% hold on
% contour(ADCP_final.time,Z,ADCP_final.u,[0,0],'Linewidth',2,'Color','k','ShowText','on','LabelSpacing',144*3);
% contour(ADCP_final.time,Z,ADCP_final.u,[0.1:0.1:1],'Color','k','ShowText','on','LabelSpacing',144*3);
% contour(ADCP_final.time,Z,ADCP_final.u,[-1:0.1:-0.1],'--k','ShowText','on','LabelSpacing',144*3);
% hold off
% ylabel('Profondeur [m]')
% title({['10°W 0°N - ZONAL VELOCITY']});
% caxis([cUaxis])
% c = colorbar;
% ylabel(c, 'U [m/s]')
% axis([min(ADCP_final.u) max(ADCP_final.u) 0 350])
% colormap(gca, parU)
% set(gca,'ydir', 'reverse');
% 
% subplot(2,1,2);
% pcolor(ADCP_final.time,Z,ADCP_final.v); shading interp
% hold on
% contour(ADCP_final.time,Z,ADCP_final.v,[0,0],'Linewidth',2,'Color','k','ShowText','on','LabelSpacing',144*3);
% contour(ADCP_final.time,Z,ADCP_final.v,[0.1:0.1:1],'Color','k','ShowText','on','LabelSpacing',144*3);
% contour(ADCP_final.time,Z,ADCP_final.v,[-1:0.1:-0.1],'--k','ShowText','on','LabelSpacing',144*3);
% hold off
% gregtick
% ylabel('Profondeur [m]')
% axis([min(ADCP_final.v) max(ADCP_final.v) 0 350])
% caxis([cVaxis])
% colormap(gca, parV)
% set(gca,'ydir', 'reverse');

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
ncid=netcdf.create([fpath,'ADCP_10W0N_2001_2019_1d.nc'],'NC_WRITE');
 
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
