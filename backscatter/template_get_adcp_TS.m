%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% template_get_adcp_TS.m
% -------------------------------
% Author : Jérémie HABASQUE - IRD
% -------------------------------
% INPUTS:
% - ADCP data processed
% - CTD file
% OUTPUTS:
% - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all

% path
addpath('..\moored_adcp_proc');

% Directory for outputs
fpath_output = '..\data_example\';

% load ADCP data
load('..\data_example\_15258_instr_01_int_filt_sub.mat');

Z = data.Z;
z = data.z_bins;
time = data.time;
freq = raw.config.sysconfig.frequency;
blen = raw.config.cell;    % bin length
nbin = raw.config.ncells;  % number of bins
bin  = 1:nbin;

% prepare S, T and sspd on ADCP-grid [d x t]
% S:   salinity on ADCP-grid;             [d x t]
% T:   in-situ temperature on ADCP-grid;  [d x t]
% sspd: soundspeed on ADCP-grid;          [d x t]

EA = data.ea(1:length(Z),:);
EA0 = 18;
tempx = data.temp;
sspdx = data.sspd;
[sspd,dummy] = meshgrid(sspdx,1:length(Z));
binlength = adcp.config.cell;
blank = adcp.config.blank;
beamangle = adcp.config.sysconfig.angle;
% mean salinity and temperature profile at deployment and recovery
file_ctd = 'Z:\PIRATA-FR25\data-processing\CTD\OS_PIRATA-FR25_CTD.nc';
lat = ncread(file_ctd, 'LATITUDE') ;
lon = ncread(file_ctd, 'LONGITUDE') ;
depth = ncread(file_ctd, 'PRES');
temp = ncread(file_ctd, 'TEMP');
sal = ncread(file_ctd, 'PSAL');
TIME = ncread(file_ctd, 'TIME') ;
TIME_CTD = TIME * 3600 *24;
for k=1:length(TIME_CTD)
    TIME_CTD(k)=datenum([2012 1 1 00 00 TIME_CTD(k)]);
end
for k=1:length(time)
    [yr, mn, dy, hr]= gregorian(time(k));
    time_adcp(k)=datenum(yr, mn, dy, hr, 00, 00);
end 
% look for the closest CTD profile
ind_profile = find(TIME_CTD>time_adcp(1), 1, 'first' );
sal_adcp = ones(length(Z),1);
temp_adcp = ones(length(Z),1);
% for each ADCP bin, find salinity and temperature values
for i_bin=1:length(Z)    
    ind_depth = find(depth(:,ind_profile)>Z(i_bin), 1, 'first' );
    sal_adcp(i_bin) = sal(ind_depth,ind_profile);
    temp_adcp(i_bin) = temp(ind_depth,ind_profile);
end
S = sal_adcp*ones(size(tempx,1),1)'; 
T = temp_adcp*ones(size(tempx,1),1)'; 

% Compute TS
out = target_strength(EA,EA0,S,T,sspd,tempx',sspdx',binlength,blank,beamangle,freq);

save([fpath_output,mooring.name, '_TS'],'out');

%% Interpolate data on a regular vertical grid
Z = fliplr(blen/2:blen:max(z(:))+blen);
Zmax = max(Z);
TS_interp = NaN(length(time),length(Z));

for i=1:length(time)
    % indice correspondant sur la grille finale Z
    ind = round((Zmax-z(1,i))/blen)+1;
    % filling the grid
    npts = min([length(Z)-ind+1 length(bin)]);
    TS_interp(i,ind:ind+npts-1) = out.ts(1:npts,i);
end
    
% Save interpolated data
save([fpath_output,mooring.name, '_TS_interp'],'out');
 
%% Figure
hf=figure('position', [0, 0, 1400, 1000]);
[C,h] = contourf(time,Z,TS_interp',50); 
set(h,'LineColor','none');
h=colorbar;
ylabel(h,'TS [dB]');
set(gca,'ydir', 'reverse');
ylabel('Depth (m)');
ylim([0,adcp.instr_depth]);
%change figure label in HH:MM
gregtick
title({[mooring.name, ' - Relative backscatter - RDI ',num2str(freq),' kHz']});

graph_name = [fpath_output, mooring.name '_relative_backscatter'];
set(hf,'Units','Inches');
pos = get(hf,'Position');
set(hf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hf,graph_name,'-dpdf','-r300');

% remove path
rmpath('.\moored_adcp_proc');
