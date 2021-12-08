%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Merge subsurface ADCP mooring data                                      %
% Autor: P. Rousselot / Date: 03/2020                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;

%% Variables
FileToMerge1     = '/media/irdcampagnes/PIRATA/PIRATA-DATA/MOORING-PIRATA-ALL/0W/merged_data/ADCP_0W0N_2016_2020_1d.nc'; %.nc file to merge with
FileToMerge2     = '/home/proussel/Documents/OUTILS/ADCP/ADCP_mooring_data_processing/FR31/0-0/ADCP_0W0N_2020_2021_1d.nc'; %.nc file to merge
step_subsampling = 1; % 1=daily
plot_data        = 1;
mooring.name     = '0W0N';
mooring.lat      = 0;
mooring.lon      = 0;
freq             = '150 khz Quartermaster';
% NCFILE info
d_fillvalue      = -9999;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read first .nc file
ncfile1.time  = ncread(FileToMerge1,'TIME');
ncfile1.depth = ncread(FileToMerge1,'DEPTH');
ncfile1.u     = ncread(FileToMerge1,'UCUR');
ncfile1.v     = ncread(FileToMerge1,'VCUR');
data.lat      = ncread(FileToMerge1,'LATITUDE');
data.lon      = ncread(FileToMerge1,'LONGITUDE');

%% Read second .nc file
ncfile2.time  = ncread(FileToMerge2,'TIME');
ncfile2.depth = ncread(FileToMerge2,'DEPTH');
ncfile2.u     = ncread(FileToMerge2,'UCUR');
ncfile2.v     = ncread(FileToMerge2,'VCUR');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create depth matrix
%depth = max(vertcat(ncfile1.depth,ncfile2.depth)):ncfile1.depth(2)-ncfile1.depth(1):min(vertcat(ncfile1.depth,ncfile2.depth));
%depth = fliplr(depth);
max(vertcat(ncfile1.depth,ncfile2.depth))
min(vertcat(ncfile1.depth,ncfile2.depth))
depth = 0:5:300;

%% Create time matrix
time  = vertcat(ncfile1.time,ncfile2.time);
[YY,MM,DD,hh,mm,ss] = datevec(time+datenum(1950,01,01));
time = julian(YY,MM,DD,hh,mm,ss);

%% Create u/v matrix
for ii = 1:length(ncfile1.time)
    u1(ii,:) = interp1(ncfile1.depth,ncfile1.u(ii,:),depth);
    v1(ii,:) = interp1(ncfile1.depth,ncfile1.v(ii,:),depth);
end
for ii = 1:length(ncfile2.time)
    u2(ii,:) = interp1(ncfile2.depth,ncfile2.u(ii,:),depth);
    v2(ii,:) = interp1(ncfile2.depth,ncfile2.v(ii,:),depth);
end

u     = vertcat(u1,u2);
v     = vertcat(v1,v2);

time  = time*ones(1,length(depth));


%% Plot data
if plot_data
    niv_u = (-1.5:0.05:0.9);
    niv_v = (-1.8:0.05:1.5);
    depth = depth.*ones(length(time),1);
    date1 = datestr(ncfile1.time(1)+datenum(1950,01,01),'YYYY');
    date2 = datestr(ncfile2.time(end)+datenum(1950,01,01),'YYYY');
    
    hf=figure('position', [0, 0, 1400, 1000]);
    %u
    subplot(2,1,1);
    colormap jet
    [C,h] = contourf(time,depth,u,niv_u);
    set(h,'LineColor','none');
    caxis(niv_u([1 end]));
    h=colorbar;
    ylabel(h,'U [m s^-^1]');
    set(gca,'ydir', 'reverse');
    ylabel('Depth (m)');
    ylim([0,max(max(depth))]);
    gregtick;
    title({[mooring.name, ' - ' num2str(date1) ' to ' num2str(date2) ' - ZONAL VELOCITY - RDI ',num2str(freq),' kHz (filtered from tide)']});
    
    %v
    subplot(2,1,2);
    [C,h] = contourf(time,depth,v,niv_v);
    set(h,'LineColor','none');
    caxis(niv_v([1 end]));
    h     = colorbar;
    ylabel(h,'V [m s^-^1]');
    set(gca,'ydir', 'reverse');
    ylabel('Depth (m)');
    ylim([0,max(max(depth))]);
    gregtick;
    title({[mooring.name, ' - ' num2str(date1) ' to ' num2str(date2) ' - MERIDIONAL VELOCITY - RDI ',num2str(freq),' kHz (filtered from tide)']});
    
    graph_name = ['ADCP_' mooring.name '_' num2str(date1) '_' num2str(date2) '_U_V_daily'];
    set(hf,'Units','Inches');
    pos = get(hf,'Position');
    set(hf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hf,graph_name,'-dpdf','-r300');
end


%%  Write netcdf file
disp('****')
disp('Creating .nc file')

% Input parameters for NETCDF Global Attributes
tc_globAttFilename      = fullfile('tools/input_GlobalAttrParameters.xls'); % JLL 2020/12 Il serait judicieux de remonter cette valeur en dÃ©but du script template_get_adcp_data.m

%% Prepare informations and variables required to create NETCDF file %% 
time       = time(:,1);
[yr_start , ~, ~]       = gregorian(time(1));
[yr_end,  ~, ~]         = gregorian(time(end));

% Read inputs metadata required for NETCDF Global Attributes
[~,~,cell_ncAttributes] = xlsread(tc_globAttFilename);

% Complete output path and filename 
tc_ncFilenam_out        = fullfile(['ADCP_',mooring.name,'_',num2str(yr_start),'_',num2str(yr_end),'_1d.nc']);

% Assign a "4D-size" (TIME,DEPTH,LATITUDE,LONGITUDE) to the ADCP current variables : UINTFILT, VINTFILT
td_uADCP                = ones(numel(time),numel(depth(1,:)),numel(data.lat),numel(data.lon)) * d_fillvalue;
td_uADCP(:,:,1,1)       = u;
td_vADCP                = ones(numel(time),numel(depth(1,:)),numel(data.lat),numel(data.lon)) * d_fillvalue;
td_vADCP(:,:,1,1)       = v;

% Flip for convention
data.Z                  = depth(1,:);
%td_uADCP                = fliplr(td_uADCP);
%td_vADCP                = fliplr(td_vADCP);

% Group general ADCP mooring informations and ADCP data to be written in NETCDF file format
struct_dataADCP         = struct('mooringName', mooring.name, 'mooringLat', mooring.lat,...
    'mooringLon', mooring.lon, 'time', time, 'depth', data.Z,...
    'u', td_uADCP, 'v', td_vADCP);

%%  Write netcdf file  %%        
f_w_ADCP_ncOS(tc_ncFilenam_out,cell_ncAttributes,struct_dataADCP,d_fillvalue);
disp('****')