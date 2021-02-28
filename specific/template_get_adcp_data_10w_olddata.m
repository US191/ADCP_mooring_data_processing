%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% template_get_adcp_data.m
% -------------------------------
% Author : Pierre ROUSSELOT - IRD (pierre.rousselot@ird.fr)
%          Jeremie HABASQUE - IRD (jeremie.habasque@ird.fr)  
% -------------------------------
% INPUTS:
% - binary raw file with .000 extension
% OUTPUTS:
% - U and V fields interpolated on a regulard grid, filtered and subsampled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear
addpath(genpath('../ADCP_mooring_data_processing'));
addpath('/home/proussel/Documents/OUTILS/TOOLS/nansuite'); % NaNSuitePath

%% META information:
% Location rawfile
rawfile          = 'FR10W000.000';        % binary file with .000 extension
fpath_output     = './10W/';                                                             % Output directory

% Cruise/mooring info
cruise.name      = 'PIRATA';                                                         % cruise name
mooring.name     = '10W0N';                                                                % '0N10W'
mooring.lat      = '00°00.000';                                                           % latitude [°']
mooring.lon      = '-10°00.000';                                                           % longitude [°']
clock_drift      = 0;                                                                   % [seconds]

% ADCP info
adcp.sn          = 500;                                                                  % ADCP serial number
adcp.type        = '300 khz';                                                             % Type : Quartermaster, longranger
adcp.direction   = 'up';                                                                  % upward-looking 'up', downward-looking 'dn'
adcp.instr_depth = 100;                                                                   % nominal instrument depth
instr            = 1;                                                                     % this is just for name convention and sorting of all mooring instruments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Convert variables
latDegInd               = strfind(mooring.lat,'°');
lonDegInd               = strfind(mooring.lon,'°');
mooring.lat             = str2double(mooring.lat(1:latDegInd-1))+str2double(mooring.lat(latDegInd+1:end-1))/60;
mooring.lon             = str2double(mooring.lon(1:lonDegInd-1))+str2double(mooring.lon(lonDegInd+1:end-1))/60;
clock_drift             = clock_drift/3600;  % convert into hrs

%% Read rawfile
disp('****')
raw_file                = [fpath_output, mooring.name '_' num2str(adcp.sn) '_instr_' sprintf('%02d',instr) '_raw.mat'];
% if exist(raw_file)
%     fprintf('Read %s\n', raw_file);
%     load(raw_file)
% else
%     fprintf('Read %s\n', rawfile);
%     raw                 = read_os3(rawfile,'all');
%     save(raw_file,'raw','-v7.3');
% end
load('UV_hour_10w_int10.mat')
d0=datenum(2004,2,6,0,0,0);
d1=datenum(2005,6,17,22,0,0); 
adcp_2003_time=d0:1/12:d1;
adcp_2003_time = adcp_2003_time';
u2 = uvmat(:,1:5976)/100;
v2 = uvmat(:,5977:5976*2)/100;
adcp_2003_bin_length = 4;
dpt1 = 4:adcp_2003_bin_length:104;
[YY,MM,DD,hh,mm,ss] = datevec(adcp_2003_time);
time=julian(YY,MM,DD,hh,mm,ss);




% Data structure
u_interp              = u2';
v_interp              = v2';
ts_interp = NaN;
data.time           = time;
Z                   = dpt1;
reply_ts = 0;

%% Horizontal interpolation, filtering and subsampling
disp('****')
horz_interp        = 1;
filt        = 1;
sub        = 1;


[uintfilt,vintfilt,tsintfilt,inttim,utid_baro,vtid_baro] = adcp_filt_sub(data,u_interp',v_interp',ts_interp',1:length(Z),reply_ts, filt, sub);
% saveas(figure(5),[fpath_output,mooring.name,'_',num2str(adcp.sn),'_instr_',num2str(instr),'_','data_raw_filt_subsampled_1'],'fig')
% saveas(figure(6),[fpath_output,mooring.name,'_',num2str(adcp.sn),'_instr_',num2str(instr),'_','data_raw_filt_subsampled_2'],'fig')

if horz_interp == 0    
    uintfilt      = u_interp';
    vintfilt      = v_interp';
    inttim        = data.time;
end

%% Save interpolated data
disp('****')
bin_start           = input('Determine first bin indice with good interpolated data: '); % bin indice where good interpolated data for the whole dataset start
disp('****')
bin_end             = input('Determine last bin indice with good interpolated data: '); % bin indice where good interpolated data for the whole dataset start
data.uintfilt       = uintfilt(bin_start:bin_end,:);
data.vintfilt       = vintfilt(bin_start:bin_end,:);
if reply_ts == 1
    data.tsintfilt  = tsintfilt(bin_start:bin_end,:);
end
data.Z              = Z(bin_start:bin_end);
data.inttim         = inttim;
%save([fpath_output, mooring.name '_' num2str(adcp.sn) '_instr_' sprintf('%02d',instr) '_int_filt_sub.mat'],'adcp','mooring','data');

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
if filt
    title({['10W0N - ZONAL VELOCITY - RDI 150 kHz (filtered from tide)']});
else
    title({[mooring.name, ' - ZONAL VELOCITY - RDI ',num2str(freq),' kHz']});
end

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
if filt
    title({['10W0N - MERIDIONAL VELOCITY - RDI 150 kHz (filtered from tide)']});
else
    title({[mooring.name, ' - MERIDIONAL VELOCITY - RDI ',num2str(freq),' kHz']});
end

graph_name = ['10W0N_nb150khz_U_V_int_filt_sub'];
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
    ylim([0,round(max(Z))]);
    %change figure label in HH:MM
    gregtick;
    if filt
        title({[mooring.name, ' - TARGET STRENGTH - RDI ',num2str(freq),' kHz (filtered from tide)']});
    else
        title({[mooring.name, ' - TARGET STRENGTH - RDI ',num2str(freq),' kHz']});
    end
    
    graph_name = [fpath_output, mooring.name, '_nb150khz_TS_int_filt_sub'];
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

ncid     = netcdf.create(['ADCP_10W0N_',num2str(yr_start),'_',num2str(yr_end),'_1d.nc'],'NC_WRITE');

%create dimension
dimidt   = netcdf.defDim(ncid,'time',length(inttim));
dimidz   = netcdf.defDim(ncid,'depth',length(data.Z));
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
netcdf.putVar(ncid,depth_ID,data.Z);
%Then store my main variable
netcdf.putVar(ncid,u_ID,data.uintfilt');
netcdf.putVar(ncid,v_ID,data.vintfilt');
if reply_ts == 1
    netcdf.putVar(ncid,ts_ID,data.tsintfilt');
end
%We're done, close the netcdf
netcdf.close(ncid);
disp('****')
% -------------------------------------------------------------------------------------------
