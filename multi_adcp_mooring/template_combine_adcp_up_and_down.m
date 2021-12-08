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
fpath = '/media/irdcampagnes/PIRATA/PIRATA-FR22/data-adjusted/MOUILLAGE_ADCP';

% Directory for outputs
fpath_output = '/media/irdcampagnes/PIRATA/PIRATA-FR22/data-adjusted/MOUILLAGE_ADCP/';
d_fillvalue = -9999;
% Contour levels for u anv v fields
niv_u = (-1.5:0.1:1.5);
niv_v = (-0.5:0.1:0.5);

%% combine up and down
load([fpath, '/up/10W0N_508_instr_01_int_filt_sub.mat']);
adcp_up=adcp;
up_u=data.uintfilt;
up_v=data.vintfilt;
up_z=data.Z;
up_time=data.inttim;
%npts_up= data.npts_up;
freq_up=adcp_up.config.sysconfig.frequency;

load([fpath, '/down/10W0N_509_instr_01_int_filt_sub.mat']);
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
depth_down_ADCP = up_z(end) + ceil(distance_between_up_and_down);
%interval grid between up and down ADCP
interval_grid = max(up_z)+adcp_up.config.cell/2-1:adcp_up.config.cell/2:min(depth_down_ADCP)-adcp_up.config.cell/2+1;
%on applique l'offset sur les profondeurs de l'ADCP down
offset_ADCP_down = min(down_z) - min(depth_down_ADCP);
down_z = down_z - offset_ADCP_down;
% down_z = flip(down_z);
% down_u = flip(down_u);
% down_v = flip(down_v);

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
Z = [up_z interval_grid down_z];
inttim = down_time;
bin_start = 1;
bin_end = length(Z);
uintfilt = u_final;
vintfilt = v_final;
data.uintfilt = uintfilt;
data.vintfilt = vintfilt;
data.inttim   = inttim;
data.Z        = Z;
mooring.name  = '10W0N';
save([fpath_output, mooring.name '_instr_int_filt_sub_up_down.mat'],'mooring','data');

%% Figure
niv_u               = (-1:0.05:1);
niv_v               = (-1:0.05:1);
close all
hf=figure('position', [0, 0, 1400, 1000]);
%u
subplot(2,1,1);
colormap jet
[C,h] = contourf(inttim,Z(bin_start:bin_end),uintfilt(bin_start:bin_end,:),niv_u);
set(h,'LineColor','none');
caxis(niv_u([1 end]));
h=colorbar;
ylabel(h,'U [m s^-^1]');
set(gca,'ydir', 'reverse');
ylabel('Depth (m)');
ylim([0,round(max(Z))]);
%change figure label in HH:MM
gregtick;
title({[mooring.name, ' - ZONAL VELOCITY - RDI 2 x 300 kHz (filtered from tide)']});


%v
subplot(2,1,2);
[C,h] = contourf(inttim,Z(bin_start:bin_end),vintfilt(bin_start:bin_end,:),niv_v);
set(h,'LineColor','none');
caxis(niv_v([1 end]));
h     = colorbar;
ylabel(h,'V [m s^-^1]');
set(gca,'ydir', 'reverse');
ylabel('Depth (m)');
ylim([0,round(max(Z))]);
%change figure label in HH:MM
gregtick;
title({[mooring.name, ' - ZONAL VELOCITY - RDI 2 x 300 kHz (filtered from tide)']});


graph_name = [fpath_output, mooring.name '_up_down_U_V_int_filt_sub'];
set(hf,'Units','Inches');
pos        = get(hf,'Position');
set(hf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(hf,graph_name,'-dpdf','-r300');

%%  Write netcdf file
disp('****')
disp('Creating .nc file')

% Input parameters for NETCDF Global Attributes
tc_globAttFilename      = fullfile('tools/input_GlobalAttrParameters.xls'); % JLL 2020/12 Il serait judicieux de remonter cette valeur en dÃ©but du script template_get_adcp_data.m

%% Prepare informations and variables required to create NETCDF file %% 
[yr_start , ~, ~]       = gregorian(data.inttim(1));
[yr_end,  ~, ~]         = gregorian(data.inttim(length(data.inttim)));

% Read inputs metadata required for NETCDF Global Attributes
[~,~,cell_ncAttributes] = xlsread(tc_globAttFilename);

% Complete output path and filename 
tc_ncFilenam_out        = fullfile(fpath_output,['ADCP_',mooring.name,'_',num2str(yr_start),'_',num2str(yr_end),'_1d.nc']);

% Assign a "4D-size" (TIME,DEPTH,LATITUDE,LONGITUDE) to the ADCP current variables : UINTFILT, VINTFILT
td_uADCP                = ones(numel(data.inttim),numel(data.Z),numel(data.lat),numel(data.lon)) * d_fillvalue;
td_uADCP(:,:,1,1)       = data.uintfilt';
td_vADCP                = ones(numel(data.inttim),numel(data.Z),numel(data.lat),numel(data.lon)) * d_fillvalue;
td_vADCP(:,:,1,1)       = data.vintfilt';

% Flip for convention
data.Z                  = fliplr(data.Z);
td_uADCP                = fliplr(td_uADCP);
td_vADCP                = fliplr(td_vADCP);

% Group general ADCP mooring informations and ADCP data to be written in NETCDF file format
struct_dataADCP         = struct('mooringName', mooring.name, 'mooringLat', mooring.lat,...
    'mooringLon', mooring.lon, 'time', data.inttim, 'depth', data.Z,...
    'u', td_uADCP, 'v', td_vADCP);

%%  Write netcdf file  %%        
f_w_ADCP_ncOS(tc_ncFilenam_out,cell_ncAttributes,struct_dataADCP,d_fillvalue);
disp('****')

