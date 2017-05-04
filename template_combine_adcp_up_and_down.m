%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% template_combine_adcp_up_and_down.m
% -------------------------------
% Author : Jérémie HABASQUE - IRD
% -------------------------------
% INPUTS:
% - ADCP up data processed (output of template_get_adcp_data.m)
% - ADCP down data processed (output of template_get_adcp_data.m) 
% OUTPUTS:
% - U and V fields interpolated on a regulard grid, filtered and subsampled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;

% path
addpath('.\moored_adcp_proc');

% Location of ADCP up and down data
fpath = '.\data_example\up_and_down\';

% Directory for outputs
fpath_output = '.\data_example\up_and_down\';

% Contour levels for u anv v fields
niv_u = (-1.5:0.1:1.5);
niv_v = (-0.5:0.1:0.5);

%% combine up and down
load([fpath, 'FR22-10W0N_508_instr_01_int_filt_sub']);
adcp_up=adcp;
up_u=data.uintfilt;
up_v=data.vintfilt;
up_z=data.Z;
up_time=data.inttim;
npts_up= data.npts_up;
freq_up=adcp_up.config.sysconfig.frequency;

load([fpath, 'FR22-10W0N_509_instr_02_int_filt_sub_down']);
adcp_down=adcp;
down_u=data.uintfilt;
down_v=data.vintfilt;
down_z=data.Z;
down_time=data.inttim;
freq_down = adcp_down.config.sysconfig.frequency;

% distance between deepest measurement of the upward looking and 
% shallowest measurement of the downward looking ADCP 
% half distance of 1st bin in upward ADCP + 
% blank of upward ADCP + 
% distance between instruments + 
% blank of downward ADCP +
% half distance of 1st bin in downward ADCP.
distance_between_instruments = 3;
distance_between_up_and_down = adcp_up.config.cell/2 + adcp_up.config.blank + distance_between_instruments + adcp_down.config.blank + adcp_down.config.cell/2;
%calculate real minimum down ADCP depth
up_z = fliplr(up_z);
depth_down_ADCP = up_z(min(npts_up)+1) + ceil(distance_between_up_and_down);
%interval grid between up and down ADCP
interval_grid = max(up_z)+adcp_up.config.cell/2-1:adcp_up.config.cell/2:min(depth_down_ADCP)-adcp_up.config.cell/2+1;
%on applique l'offset sur les profondeurs de l'ADCP down
offset_ADCP_down = min(down_z) - min(depth_down_ADCP);
down_z = down_z - offset_ADCP_down;

% initialisation de la matrice combinée up + down
[B,a] = size(down_u); [B1,a1] = size(up_u);
BB = B+B1+1 ;
DOWN = NaN(BB,a);
V_DOWN = NaN(BB,a);
u_final = NaN(BB,a);
v_final = NaN(BB,a);

% on alimente la matrice avec les données up
u_final(1:B1,1:a1) = up_u(B1:-1:1,:);
v_final(1:B1,1:a1) = up_v(B1:-1:1,:);

for i = 1:length(down_u)
    u_final(B1+length(interval_grid)+1:B1+length(interval_grid)+B,i) = down_u(:,i);
    v_final(B1+length(interval_grid)+1:B1+length(interval_grid)+B,i) = down_v(:,i);
end

%% figure finale
z_final = [up_z interval_grid down_z];

hf=figure('position', [0, 0, 1400, 1000]);
%u
subplot(2,1,1);
[~,h] = contourf(down_time,z_final,u_final,niv_u); 
set(h,'LineColor','none');
caxis(niv_u([1 end]));
h=colorbar;
ylabel(h,'U [m s^-^1]');
set(gca,'ydir', 'reverse');
ylabel('Depth (m)');
gregtick
title({[mooring.name, ' - ZONAL VELOCITY - ADCP UP : RDI ',num2str(freq_up),' kHz - ADCP DOWN : RDI ',num2str(freq_down),' kHz']});

%v
subplot(2,1,2);
[C,h] = contourf(down_time,z_final,v_final,niv_v); 
set(h,'LineColor','none');
caxis(niv_v([1 end]));
h=colorbar;
ylabel(h,'V [m s^-^1]');
set(gca,'ydir', 'reverse');
ylabel('Depth (m)');
gregtick
title({[mooring.name, ' - MERIDIONAL VELOCITY - ADCP UP : RDI ',num2str(freq_up),' kHz - ADCP DOWN : RDI ',num2str(freq_down),' kHz']});

% Save combined data
data.time = down_time;
data.u_final = u_final;
data.v_final = v_final;
data.z_final = z_final;
save([fpath, mooring.name '_UP_DOWN_int_filt_sub.mat'],'adcp','mooring','data','raw');
 
graph_name = [fpath_output, mooring.name '_U_V_UP_DOWN_int_filt_sub'];
set(hf,'Units','Inches');
pos = get(hf,'Position');
set(hf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hf,graph_name,'-dpdf','-r300');

% rmpath
rmpath('.\moored_adcp_proc');
clear all; close all;

%%  Write netcdf file           
[yr_start , ~, ~] = gregorian(inttim(1));
[yr_end,  ~, ~] = gregorian(inttim(length(inttim)));

ncid=netcdf.create([fpath_output,'ADCP_',mooring.name, '_',num2str(yr_start),'_',num2str(yr_end),'_1d.nc'],'NC_WRITE');
 
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
netcdf.putVar(ncid,time_ID,down_time);
netcdf.putVar(ncid,depth_ID,z_final);  
%Then store my main variable
netcdf.putVar(ncid,u_ID,u_final);
netcdf.putVar(ncid,v_ID,v_final);
%We're done, close the netcdf
netcdf.close(ncid);

