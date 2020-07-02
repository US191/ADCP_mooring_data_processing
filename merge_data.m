%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Merge subsurface ADCP mooring data                                      %
% Autor: P. Rousselot / Date: 03/2020                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;

%% Variables
FileToMerge1     = 'C:\Users\proussel\Documents\OUTILS\ADCP\ADCP_mooring_data_processing\FR28\ADCP_0N0W_2016_2018_1d.nc'; %.nc file to merge with
FileToMerge2     = 'C:\Users\proussel\Documents\OUTILS\ADCP\ADCP_mooring_data_processing\FR30\ADCP_0N0W_2018_2020_1d.nc'; %.nc file to merge
FileToMerge1     = '/home/proussel/Documents/OUTILS/ADCP/ADCP_mooring_data_processing/IDMX/ADCP_0N130E_2013_2013_1d.nc'; %.nc file to merge with
FileToMerge2     = '/home/proussel/Documents/OUTILS/ADCP/ADCP_mooring_data_processing/IDMX/ADCP_0N130E_2013_2016_1d.nc'; %.nc file to merge
step_subsampling = 1; % 1=daily
plot_data        = 0;
mooring.name     = '0N0E';
freq             = '150';
mooring.name     = '0N130E';
freq             = '75';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read first .nc file
ncfile1.time  = ncread(FileToMerge1,'time');
ncfile1.depth = ncread(FileToMerge1,'depth');
ncfile1.u     = ncread(FileToMerge1,'u');
ncfile1.v     = ncread(FileToMerge1,'v');

%% Read second .nc file
ncfile2.time  = ncread(FileToMerge2,'time');
ncfile2.depth = ncread(FileToMerge2,'depth');
ncfile2.u     = ncread(FileToMerge2,'u');
ncfile2.v     = ncread(FileToMerge2,'v');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create depth matrix
depth = max(vertcat(ncfile1.depth,ncfile2.depth)):ncfile1.depth(2)-ncfile1.depth(1):min(vertcat(ncfile1.depth,ncfile2.depth));
depth = fliplr(depth);

%% Create time matrix
time  = vertcat(ncfile1.time,ncfile2.time);
time  = time*ones(1,length(depth));

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

%% Plot data
if plot_data
    niv_u = (-1.5:0.05:0.9);
    niv_v = (-1.8:0.05:1.5);
    depth = depth.*ones(length(time),1);
    date1 = julian(ncfile1.time(1));
    date2 = julian(ncfile2.time(end));
    
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
    title({[mooring.name, ' - ' num2str(date1(1)) ' to ' num2str(date2(1)) ' - ZONAL VELOCITY - RDI ',num2str(freq),' kHz (filtered from tide)']});
    
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
    title({[mooring.name, ' - ' num2str(date1(1)) ' to ' num2str(date2(1)) ' - MERIDIONAL VELOCITY - RDI ',num2str(freq),' kHz (filtered from tide)']});
    
    graph_name = ['ADCP_' mooring.name '_' num2str(date1(1)) '_' num2str(date2(1)) '_U_V_daily'];
    set(hf,'Units','Inches');
    pos = get(hf,'Position');
    set(hf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hf,graph_name,'-dpdf','-r300');
end

%%  Write netcdf file
time       = time(:,1);
timed      = gregorian(time);
timed(:,5) = 0;
timed(:,6) = 0;
time       = datenum(timed)-datenum(1950,01,01);

disp('****')
disp('Creating .nc file')
yr_start = datestr(time(1)+datenum(1950,01,01), 'yyyy');
yr_end   = datestr(time(end)+datenum(1950,01,01), 'yyyy');

ncid     = netcdf.create(['ADCP_',mooring.name,'_',num2str(yr_start),'_',num2str(yr_end),'.nc'],'NC_CLOBBER');

%create dimension
dimidz   = netcdf.defDim(ncid, 'DEPTH', length(depth));
dimidt   = netcdf.defDim(ncid, 'TIME', netcdf.getConstant('NC_UNLIMITED'));
%Define IDs for the dimension variables (pressure,time,latitude,...)
time_ID  = netcdf.defVar(ncid,'TIME','double',dimidt);
netcdf.putAtt(ncid, time_ID, 'long_name', 'Time');
netcdf.putAtt(ncid, time_ID, 'axis', 'T');
netcdf.putAtt(ncid, time_ID, 'units', 'days since 1950-01-01T00:00:00Z');
netcdf.putAtt(ncid, time_ID, 'time_origin', '1-JAN-1950:00:00:00');
depth_ID = netcdf.defVar(ncid,'DEPTH','double',dimidz);
netcdf.putAtt(ncid, depth_ID, 'long_name', 'Depth');
netcdf.putAtt(ncid, depth_ID, 'axis', 'Z');
netcdf.putAtt(ncid, depth_ID, 'units', 'meters');
netcdf.putAtt(ncid, depth_ID, 'positive', 'down');
netcdf.putAtt(ncid, depth_ID, 'point_spacing', 'even');
%Define the main variable ()
u_ID     = netcdf.defVar(ncid,'U','double',[dimidz dimidt]);
netcdf.putAtt(ncid, u_ID, 'long_name', 'Zonal velocities');
v_ID     = netcdf.defVar(ncid,'V','double',[dimidz dimidt]);
netcdf.putAtt(ncid, v_ID, 'long_name', 'Meridional velocities');
%We are done defining the NetCdf
netcdf.endDef(ncid);
%Then store the dimension variables in
netcdf.putVar(ncid, depth_ID, depth(1,:));
netcdf.putVar(ncid, time_ID, 0, length(time), time);
%Then store my main variable
netcdf.putVar(ncid, u_ID, u');
netcdf.putVar(ncid, v_ID, v');
%We're done, close the netcdf
netcdf.close(ncid);
disp('****')