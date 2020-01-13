load('FR28/0N0W_15258_instr_05_int_filt_sub.mat','adcp','mooring','data');
addpath(genpath('../ADCP_mooring_data_processing'));
addpath('/home/proussel/Documents/OUTILS/TOOLS/nansuite'); % NaNSuitePath
fpath_output='FR28/';
reply_ts=1;
freq=150;
Z = data.Z;
inttim = data.inttim;
bin_start = 1;
bin_end = length(Z);
uintfilt = data.uintfilt;
vintfilt = data.vintfilt;
tsintfilt = data.tsintfilt;

%% Figure
niv_u               = (-1:0.05:1);
niv_v               = (-1:0.05:1);
close all
hf=figure('position', [0, 0, 1400, 1500]);
%u
subplot(3,1,1);
colormap jet
[C,h] = contourf(inttim,Z(bin_start:bin_end),uintfilt(bin_start:bin_end,:),niv_u); 
set(h,'LineColor','none');
caxis(niv_u([1 end]));
h=colorbar;
ylabel(h,'U [m s^-^1]');
set(gca,'ydir', 'reverse');
ylabel('Depth (m)');
ylim([0,adcp.instr_depth]);
%change figure label in HH:MM
gregtick;
title({[mooring.name, ' - ZONAL VELOCITY - RDI ',num2str(freq),' kHz']});

%v
subplot(3,1,2);
[C,h] = contourf(inttim,Z(bin_start:bin_end),vintfilt(bin_start:bin_end,:),niv_v); set(h,'LineColor','none');
caxis(niv_v([1 end]));
h     = colorbar;
ylabel(h,'V [m s^-^1]');
set(gca,'ydir', 'reverse');
ylabel('Depth (m)');
ylim([0,adcp.instr_depth]);
%change figure label in HH:MM
gregtick;
title({[mooring.name, ' - MERIDIONAL VELOCITY - RDI ',num2str(freq),' kHz']});


%wind
subplot(3,1,3);
wind = dlmread('/home/proussel/Documents/OUTILS/ADCP/ADCP_mooring_data_processing/FR28/w0n0e_dy.ascii');
x = (1:5:684)';
y = zeros(size(wind(1:5:end,3)));
u = wind(1:5:end,3);
v = wind(1:5:end,4);
time = strcat(num2str(wind(1:5:end,1)),num2str(wind(1:5:end,2)));
for ii = 1 :length(time)
    ti(ii) = julian(str2double(time(ii,1:4)),str2double(time(ii,5:6)),str2double((time(ii,7:8))));
end
quiver(ti',y,u,v)
axis([ti(1) ti(end) 0 50])
colorbar
gregtick;
title({[mooring.name, ' - WIND']});

graph_name = [fpath_output, mooring.name '_U_V_int_filt_sub'];
set(hf,'Units','Inches');
pos        = get(hf,'Position');
set(hf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(hf,graph_name,'-dpdf','-r300');

if reply_ts == 1
hf=figure('position', [0, 0, 1400, 500]);
colormap jet
pcolor(inttim,Z(bin_start:bin_end),tsintfilt(bin_start:bin_end,:)); shading interp; 
h     = colorbar;
ylabel(h,'Target Strength [dB1]');
set(gca,'ydir', 'reverse');
ylabel('Depth (m)');
ylim([0,adcp.instr_depth]);
%change figure label in HH:MM
gregtick;
title({[mooring.name, ' - TARGET STRENGTH - RDI ',num2str(freq),' kHz']});   

graph_name = [fpath_output, mooring.name '_TS_int_filt_sub'];
set(hf,'Units','Inches');
pos        = get(hf,'Position');
set(hf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
saveas(hf,graph_name,'jpeg');
end



% %% Plot tide
% hf=figure('position', [0, 0, 1400, 1000]);
% niv_tide = (-5:0.2:5);
% utid_baro = utid_baro * 100;
% vtid_baro = vtid_baro * 100;
% %u
% subplot(2,1,1);
% colormap jet
% [C,h] = contourf(inttim,Z(bin_start:bin_end),utid_baro(bin_start:bin_end,:),niv_tide); 
% set(h,'LineColor','none');
% caxis(niv_tide([1 end]));
% h=colorbar;
% ylabel(h,'U [cm s^-^1]');
% set(gca,'ydir', 'reverse');
% ylabel('Depth (m)');
% ylim([0,adcp.instr_depth]);
% %change figure label in HH:MM
% gregtick;
% title({[mooring.name, ' - ZONAL TIDE VELOCITY - RDI ',num2str(freq),' kHz']});
% 
% %v
% subplot(2,1,2);
% [C,h] = contourf(inttim,Z(bin_start:bin_end),vtid_baro(bin_start:bin_end,:),niv_tide); 
% set(h,'LineColor','none');
% caxis(niv_tide([1 end]));
% h     = colorbar;
% ylabel(h,'V [cm s^-^1]');
% set(gca,'ydir', 'reverse');
% ylabel('Depth (m)');
% ylim([0,adcp.instr_depth]);
% %change figure label in HH:MM
% gregtick;
% title({[mooring.name, ' - MERIDIONAL TIDE VELOCITY - RDI ',num2str(freq),' kHz']});

%%  Write netcdf file     
disp('****')
disp('Creating .nc file')
[yr_start , ~, ~] = gregorian(inttim(1));
[yr_end,  ~, ~]   = gregorian(inttim(length(inttim)));

ncid     = netcdf.create([fpath_output,'ADCP_',mooring.name,'_',num2str(yr_start),'_',num2str(yr_end),'_1d.nc'],'NC_WRITE');
 
%create dimension
dimidt   = netcdf.defDim(ncid,'time',length(inttim));
dimidz   = netcdf.defDim(ncid,'depth',length(Z));
%Define IDs for the dimension variables (pressure,time,latitude,...)
time_ID  = netcdf.defVar(ncid,'time','double',dimidt);
depth_ID = netcdf.defVar(ncid,'depth','double',dimidz);
%Define the main variable ()
u_ID     = netcdf.defVar(ncid,'u','double',[dimidt dimidz]);
v_ID     = netcdf.defVar(ncid,'v','double',[dimidt dimidz]);
if reply_ts == 1
    ts_ID     = netcdf.defVar(ncid,'ts','double',[dimidt dimidz]);
end
%We are done defining the NetCdf
netcdf.endDef(ncid);
%Then store the dimension variables in
netcdf.putVar(ncid,time_ID,inttim);
netcdf.putVar(ncid,depth_ID,Z);  
%Then store my main variable
netcdf.putVar(ncid,u_ID,uintfilt');
netcdf.putVar(ncid,v_ID,vintfilt');
if reply_ts == 1
    netcdf.putVar(ncid,ts_ID,tsintfilt');
end
%We're done, close the netcdf
netcdf.close(ncid);
disp('****')
% -------------------------------------------------------------------------------------------