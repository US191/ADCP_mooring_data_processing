%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Merge subsurface ADCP mooring data                                      %
% Autor: P. Rousselot / Date: 03/2020                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;

%% Variables
FileToMerge1     = 'F:\Encours\PIRATA_ADCP-MOORINGS\Convert_ADCP/ADCP_23W0N_2001_2015.nc'; %.nc file to merge with
FileToMerge2     = 'F:\Encours\PIRATA_ADCP-MOORINGS\23W\0_Final_Files/ADCP_0N23W_2015_2016_1d_up_down.nc'; %.nc file to merge
% FileToMerge3     = '/media/irdcampagnes/PIRATA/PIRATA_ADCP-MOORINGS/10W/0_Final_Files/ADCP_10W0N_2004_2005_1d.nc'; %.nc file to merge
% FileToMerge1     = '/home/proussel/Documents/OUTILS/ADCP/ADCP_mooring_data_processing/IDMX/ADCP_0N130E_2013_2013_1d.nc'; %.nc file to merge with
% FileToMerge2     = '/home/proussel/Documents/OUTILS/ADCP/ADCP_mooring_data_processing/IDMX/ADCP_0N130E_2013_2016_1d.nc'; %.nc file to merge
step_subsampling = 1; % 1=daily
plot_data        = 1;
mooring.name     = '23W0N';
freq             = 'Variable';
% mooring.name     = '0N130E';
% freq             = '75';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read first .nc file
ncfile1.time  = ncread(FileToMerge1,'TIME');
ncfile1.depth = ncread(FileToMerge1,'DEPTH');
ncfile1.u     = ncread(FileToMerge1,'U')';
ncfile1.v     = ncread(FileToMerge1,'V')';

%% Read second .nc file
ncfile2.time  = ncread(FileToMerge2,'TIME');
ncfile2.depth = ncread(FileToMerge2,'DEPTH');
ncfile2.u     = ncread(FileToMerge2,'UCUR');
ncfile2.v     = ncread(FileToMerge2,'VCUR');

% %% Read third .nc file
% ncfile3.time  = ncread(FileToMerge3,'TIME');
% ncfile3.depth = ncread(FileToMerge3,'DEPTH');
% ncfile3.u     = ncread(FileToMerge3,'UCUR');
% ncfile3.v     = ncread(FileToMerge3,'VCUR');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create depth matrix
depth = max(vertcat(ncfile1.depth,ncfile2.depth)):ncfile1.depth(2)-ncfile1.depth(1):min(vertcat(ncfile1.depth,ncfile2.depth));
depth = fliplr(depth);
max(vertcat(ncfile1.depth,ncfile2.depth))
% min(vertcat(ncfile1.depth,ncfile2.depth))
depth = 0:5:1000;

%% Create time matrix
time  = vertcat(ncfile1.time,ncfile2.time);
[YY,MM,DD,hh,mm,ss] = datevec(time+datenum(1950,01,01));
time = julian(YY,MM,DD,hh,mm,ss);
%time  = time + datenum(1950,1,1,0,0,0) + datenum(2020,01,01);  % -0.5 is necessary to remove 12h00 introduced by datetime function starting at noon instead of midnight as performed by julian function previously used 
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
% for ii = 1:length(ncfile3.time)
%     u3(ii,:) = interp1(ncfile3.depth,ncfile3.u(ii,:),depth);
%     v3(ii,:) = interp1(ncfile3.depth,ncfile3.v(ii,:),depth);
% end
u     = vertcat(u1,u2);
v     = vertcat(v1,v2);

u = u';
v = v';

% create a continuous series of daily data, ranging from min(d) to max(d)
final_time = ceil(min(min(time))):1:floor(max(max(time)));
ADCP_final.u = NaN(length(depth),length(final_time));
ADCP_final.v = NaN(length(depth),length(final_time));

for i_time = 1:length(final_time)
    for j_time = 1:length(time)
        if final_time(i_time) == time(j_time,1)
            ADCP_final.u(:,i_time)=u(:,j_time);
            ADCP_final.v(:,i_time)=v(:,j_time);
        end
    end
end

time = final_time;
time  = time'*ones(1,length(depth));
time = time';
u    = ADCP_final.u;
v    = ADCP_final.v;

%% Plot data
if plot_data
    niv_u = (-1.5:0.05:1);
    niv_v = (-1.8:0.05:1.5);
    depth = depth.*ones(length(time),1);
    depth = depth';
    %date1 = julian(ncfile1.time(1));
    %date2 = julian(ncfile2.time(end));
    date1 = datestr(ncfile1.time(1)+datenum(1950,01,01),'yyyy');
    date2 = datestr(ncfile2.time(end)+datenum(1950,01,01),'yyyy');
    
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
    %datetick('x','mm-yyyy');
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
time       = time(1,:)';
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
dimidz   = netcdf.defDim(ncid, 'DEPTH', size(depth,2));
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
netcdf.putVar(ncid, depth_ID, depth);
netcdf.putVar(ncid, time_ID, 0, length(time), time);
%Then store my main variable
netcdf.putVar(ncid, u_ID, u);
netcdf.putVar(ncid, v_ID, v);
%We're done, close the netcdf
netcdf.close(ncid);
disp('****')