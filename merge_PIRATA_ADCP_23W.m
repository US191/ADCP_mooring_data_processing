%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% merge_PIRATA_ADCP_23W.m
% ------------------------------------
% Merge ADCP datasets from 2001 to 2016
% -------------------------------
% Author : Jérémie HABASQUE - IRD
% -------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;
% 
% r = [0 0 1];       %# start
% w = [.9 .9 .9];    %# middle
% b = [1 0 0];       %# end
% 
% %# colormap of size 64-by-3, ranging from red -> white -> blue
% c1 = zeros(32,3); c2 = zeros(32,3);
% for i=1:3
%     c1(:,i) = linspace(r(i), w(i), 32);
%     c2(:,i) = linspace(w(i), b(i), 32);
% end
% c = [c1(1:end-1,:);c2];

fpath = 'C:\Users\jhabasqu\Desktop\PIRATA\ADCP_mooring_data\23W-0N\merge_data\';
mooring.name='23W0N';

step_subsampling = 0.5; %subsampling on a regular 12-hour grid  

%% Read data

%2001-2002
load('C:\Users\jhabasqu\Desktop\PIRATA\ADCP_mooring_data\23W-0N\2001-2002 - kpo_0611\23W0N_kpo_0611_instr_01_int_filt_sub.mat');
adcp_2002_time = data.inttim;
adcp_2002_u = data.uintfilt;
adcp_2002_v = data.vintfilt;
adcp_2002_z = data.Z; 

%2004-2005
load('C:\Users\jhabasqu\Desktop\PIRATA\ADCP_mooring_data\23W-0N\2004-2005 - kpo_0612\23W0N_kpo_0612_instr_01_int_filt_sub.mat');
adcp_2004_time = data.inttim;
adcp_2004_u = data.uintfilt;
adcp_2004_v = data.vintfilt;
adcp_2004_z = data.Z; 

%2006-2008
load('C:\Users\jhabasqu\Desktop\PIRATA\ADCP_mooring_data\23W-0N\2006-2008 - kpo_1001\23W0N_kpo_1001_instr_01_int_filt_sub.mat');
adcp_2006_time = data.inttim;
adcp_2006_u = data.uintfilt;
adcp_2006_v = data.vintfilt;
adcp_2006_z = data.Z; 

%2008-2009
load('C:\Users\jhabasqu\Desktop\PIRATA\ADCP_mooring_data\23W-0N\2008-2009 - kpo_1023\23W0N_UP_DOWN_int_filt_sub.mat');
adcp_2008_time = data.inttim;
adcp_2008_u = data.u_final;
adcp_2008_v = data.v_final;
adcp_2008_z = data.z_final; 

% 2010-2011
load('C:\Users\jhabasqu\Desktop\PIRATA\ADCP_mooring_data\23W-0N\2009-2011 - kpo_1044\23W0N_UP_DOWN_int_filt_sub.mat');
adcp_2010_time = data.inttim;
adcp_2010_u = data.u_final;
adcp_2010_v = data.v_final;
adcp_2010_z = data.z_final; 

% 2011-2012
load('C:\Users\jhabasqu\Desktop\PIRATA\ADCP_mooring_data\23W-0N\2011-2012 - kpo_1063\23W0N_UP_DOWN_int_filt_sub.mat');
adcp_2011_time = data.inttim;
adcp_2011_u = data.u_final;
adcp_2011_v = data.v_final;
adcp_2011_z = data.z_final; 

% 2012-2014
load('C:\Users\jhabasqu\Desktop\PIRATA\ADCP_mooring_data\23W-0N\2012-2014 - kpo_1089\23W0N_UP_DOWN_int_filt_sub.mat');
adcp_2012_time = data.inttim;
adcp_2012_u = data.u_final;
adcp_2012_v = data.v_final;
adcp_2012_z = data.z_final; 

% 2014-2015
load('C:\Users\jhabasqu\Desktop\PIRATA\ADCP_mooring_data\23W-0N\2014-2015 - kpo_1125\23W0N_UP_DOWN_int_filt_sub.mat');
adcp_2014_time = data.inttim;
adcp_2014_u = data.u_final;
adcp_2014_v = data.v_final;
adcp_2014_z = data.z_final; 

% 2015-2016
load('C:\Users\jhabasqu\Desktop\PIRATA\ADCP_mooring_data\23W-0N\2015-2016 - kpo_1140\23W0N_UP_DOWN_int_filt_sub.mat');
adcp_2015_time = data.inttim;
adcp_2015_u = data.u_final;
adcp_2015_v = data.v_final;
adcp_2015_z = data.z_final; 

%% Interpolate data on a regular vertical grid
  
%Z = fliplr(min_depth:max_bin_length:max_depth);
Z = fliplr(0:5:350); % on se base sur la maille des donnees TACE etendue a 350m
Zmax = max(Z);

%interpolation for each timestep for 2002 data
u_interp_2002 = NaN(length(Z),length(adcp_2002_time));
v_interp_2002 = NaN(length(Z),length(adcp_2002_time));
for i=1:length(adcp_2002_time)    
     ind_ok = find(~isnan(adcp_2002_u(:,i)));
     u_interp_2002(:,i) = interp1(adcp_2002_z(ind_ok),adcp_2002_u(ind_ok,i),Z);
     v_interp_2002(:,i) = interp1(adcp_2002_z(ind_ok),adcp_2002_v(ind_ok,i),Z);    
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

%interpolation for each timestep for 2008 data
u_interp_2008 = NaN(length(Z),length(adcp_2008_time));
v_interp_2008 = NaN(length(Z),length(adcp_2008_time));
for i=1:length(adcp_2008_time)    
     ind_ok = find(~isnan(adcp_2008_u(:,i)));
     u_interp_2008(:,i) = interp1(adcp_2008_z(ind_ok),adcp_2008_u(ind_ok,i),Z);
     v_interp_2008(:,i) = interp1(adcp_2008_z(ind_ok),adcp_2008_v(ind_ok,i),Z);    
end
 
%interpolation for each timestep for 2010 data
u_interp_2010 = NaN(length(Z),length(adcp_2010_time));
v_interp_2010 = NaN(length(Z),length(adcp_2010_time));
for i=1:length(adcp_2010_time)    
     ind_ok = find(~isnan(adcp_2010_u(:,i)));
     u_interp_2010(:,i) = interp1(adcp_2010_z(ind_ok),adcp_2010_u(ind_ok,i),Z);
     v_interp_2010(:,i) = interp1(adcp_2010_z(ind_ok),adcp_2010_v(ind_ok,i),Z);    
end

%interpolation for each timestep for 2011 data
u_interp_2011 = NaN(length(Z),length(adcp_2011_time));
v_interp_2011 = NaN(length(Z),length(adcp_2011_time));
for i=1:length(adcp_2011_time)    
     ind_ok = find(~isnan(adcp_2011_u(:,i)));
     u_interp_2011(:,i) = interp1(adcp_2011_z(ind_ok),adcp_2011_u(ind_ok,i),Z);
     v_interp_2011(:,i) = interp1(adcp_2011_z(ind_ok),adcp_2011_v(ind_ok,i),Z);    
end

%interpolation for each timestep for 2012 data
u_interp_2012 = NaN(length(Z),length(adcp_2012_time));
v_interp_2012 = NaN(length(Z),length(adcp_2012_time));
for i=1:length(adcp_2012_time)    
     ind_ok = find(~isnan(adcp_2012_u(:,i)));
     u_interp_2012(:,i) = interp1(adcp_2012_z(ind_ok),adcp_2012_u(ind_ok,i),Z);
     v_interp_2012(:,i) = interp1(adcp_2012_z(ind_ok),adcp_2012_v(ind_ok,i),Z);    
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
  
%% combine all data
all_time = [adcp_2002_time adcp_2004_time adcp_2006_time adcp_2008_time adcp_2010_time adcp_2011_time adcp_2012_time adcp_2014_time adcp_2015_time];
all_u_interp = [u_interp_2002 u_interp_2004 u_interp_2006 u_interp_2008 u_interp_2010 u_interp_2011 u_interp_2012 u_interp_2014 u_interp_2015];
all_v_interp = [v_interp_2002 v_interp_2004 v_interp_2006 v_interp_2008 v_interp_2010 v_interp_2011 v_interp_2012 v_interp_2014 v_interp_2015];

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
save([fpath, 'ADCP_',mooring.name '_2001_2016_int_filt_sub.mat'],'ADCP_final');

%% FIGURES
niv_u = (-1.5:0.1:1.5);
niv_v = (-0.5:0.1:0.5);

hf=figure('position', [0, 0, 1400, 1000]);
%u
subplot(2,1,1);
[C,h] = contourf(ADCP_final.time,Z,ADCP_final.u,niv_u);
set(h,'LineColor','none');
caxis(niv_u([1 end]));
h=colorbar;
%colormap(c);
ylabel(h,'U [m s^-^1]');
set(gca,'ydir', 'reverse');
ylim([0 350]);
ylabel('Depth (m)');
%change figure label in HH:MM
gregtick
title({['23°W - ZONAL VELOCITY']});

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
title({['23°W - MERIDIONAL VELOCITY']});

graph_name = [fpath, 'ADCP_23W0N_2001_2016_U_V_daily'];
set(hf,'Units','Inches');
pos = get(hf,'Position');
set(hf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hf,graph_name,'-dpdf','-r300');

% Histogramme des valeurs U et V

hf=figure('position', [0, 0, 1400, 1000]);
subplot(1,2,1); hist(ADCP_final.u(:),100); xlabel('U [m s^-^1]');
subplot(1,2,2); hist(ADCP_final.v(:),100); xlabel('V [m s^-^1]');

graph_name = [fpath, 'ADCP_23W0N_2001_2016_U_V_histo'];
set(hf,'Units','Inches');
pos = get(hf,'Position');
set(hf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hf,graph_name,'-dpdf','-r300');

%%  Write netcdf file                           
ncid=netcdf.create([fpath,'ADCP_23W0N_2001_2016_1d.nc'],'NC_WRITE');
 
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

