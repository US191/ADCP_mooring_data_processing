%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% template_get_adcp_data.m
% -------------------------------
% Author : Jeremie HABASQUE / Pierre Rousselot - IRD
% -------------------------------
% INPUTS:
% - binary raw file with .000 extension
% OUTPUTS:
% - U and V fields interpolated on a regulard grid, filtered and subsampled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear
addpath(genpath('../ADCP_mooring_data_processing'));
addpath('/home/proussel/Documents/OUTILS/TOOLS/nansuite'); % NaNSuitePath

% First part --------------------------------------------------------------------------------------------------------------------
%% META information:
% Location rawfile
rawfile          = './FR26_000.000';        % binary file with .000 extension
fpath_output     = './FR28/';               % Output directory
 
% Cruise/mooring info
cruise.name      = 'PIRATA-FR28';           % cruise name
mooring.name     = '0N0W';                  % '0N10W'
mooring.lat      = '00°01.060';             % latitude [°']
mooring.lon      = '-000°00.330';           % longitude [°']
clock_drift      = 0/3600;                  % [seconds]

% ADCP info
adcp.sn          = 15258;                   % ADCP serial number
adcp.type        = '150 khz Quartermaster'; % Type : Quartermaster, longranger
adcp.direction   = 'up';                    % upward-looking 'up', downward-looking 'dn'
adcp.instr_depth = 300;                     % nominal instrument depth
instr            = 1;                       % this is just for name convention and sorting of all mooring instruments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Convert variables
latDegInd               = strfind(mooring.lat,'°');  
lonDegInd               = strfind(mooring.lon,'°');  
mooring.lat             = str2double(mooring.lat(1:latDegInd-1))+str2double(mooring.lat(latDegInd+1:end-1))/60;
mooring.lon             = str2double(mooring.lon(1:lonDegInd-1))+str2double(mooring.lon(lonDegInd+1:end-1))/60;
clock_drift             = clock_drift/3600;  % convert into hrs

%% Read rawfile
disp('****')
fprintf('Read %s\n', rawfile);
raw                     = read_os3(rawfile,'all');

%% Correct clock drift
time0                   = julian(raw.juliandate);
clockd                  = linspace(0, clock_drift, length(time0));
raw.juliandate          = raw.juliandate - clockd / 24;  
disp('****')
disp('Correct clock drift')

%% Plot pressure and temperature sensor
figure;
subplot(2,1,1)
plot(raw.pressure);
detrend_sdata = detrend(raw.pressure);
trend         = raw.pressure - detrend_sdata;
hold on
plot(trend, 'r--')
hold off
title('Pressure sensor');
ylabel('Depth [m]');
xlabel('Time index');
grid on; 
saveas(gcf,[fpath_output,mooring.name,'_',num2str(adcp.sn),'_instr_',num2str(instr),'_','Pressure_sensor'],'fig')

subplot(2,1,2)
plot(raw.temperature);
title('Temperature sensor');
ylabel('Temperature [°C]');
xlabel('Time index');
grid on; 
% Second part --------------------------------------------------------------------------------------------------------------------
%% Determine first and last indiced when instrument was at depth (you can do this by plotting 'raw.pressure' for example        
disp('****')
first               = input('Determine first indice when instrument was at depth (with pres/temp plot): ');
disp('****')
last                = input('Determine last indice when instrument was at depth (with pres/temp plot): ');

%% amplitude of the bins / Correction ADCP's depth
ea                  = squeeze(mean(raw.amp(:,:,first:last),2));   
figure; 
colormap jet; 
pcolor(ea); 
shading flat; 
title('Amplitude of the bins'); colorbar;
ylabel('Bins');
xlabel('Time index'); 
saveas(gcf,[fpath_output,mooring.name,'_',num2str(adcp.sn),'_instr_',num2str(instr),'_','Amplitude_bins'],'fig')

% Third part --------------------------------------------------------------------------------------------------------------------

%% If upward looking: determine range of surface bins used for instrument depth correction below!
disp('****')
sbins               = input('Determine range of surface bins used for instrument depth correction (with aplitude plot, ie. 30:35): ');

%% Correct velocity data with external T/S sensor

%% Extract data
freq                = raw.config.sysconfig.frequency;
u2                  = squeeze(raw.vel(:,1,first:last));
v2                  = squeeze(raw.vel(:,2,first:last));
w                   = squeeze(raw.vel(:,3,first:last));
vel_err             = squeeze(raw.vel(:,4,first:last));  % the difference in vertical velocity between the two pairs of transducers 
time                = raw.juliandate(first:last);
ang                 = [raw.pitch(first:last) raw.roll(first:last) raw.heading(first:last)]; 
soundspeed          = raw.soundspeed(first:last);
temp                = raw.temperature(first:last);
press               = raw.pressure(first:last);
if press < 0
    press = -press;
end
nbin                = raw.config.ncells;      % number of bins
bin                 = 1:nbin;                 % bin number
blen                = raw.config.cell;        % bin length
tlen                = raw.config.pulse;       % pulse length
lag                 = raw.config.lag;         % transmit lag distance
blnk                = raw.config.blank;       % blank distance after transmit
fbind               = raw.config.bin1distance;% middle of first bin distance
dt                  = (time(2)-time(1))*24;   % Sampling interval in hours

%% Calculate Magnetic deviation values
[a,~]                   = gregorian(raw.juliandate(1));
magnetic_deviation_ini  = magdev(mooring.lat,mooring.lon,0,a+(raw.juliandate(1)-julian(a,1,1,0,0,0))/365.25);
[a,~]                   = gregorian(raw.juliandate(end));
magnetic_deviation_end  = magdev(mooring.lat,mooring.lon,0,a+(raw.juliandate(end)-julian(a,1,1,0,0,0))/365.25);
rot                     = (magnetic_deviation_ini+magnetic_deviation_end)/2;  
mag_dev                 = linspace(magnetic_deviation_ini, magnetic_deviation_end, length(time0)); 
mag_dev             = mag_dev(first:last);

%% Correction of magnetic deviation
disp('****')
disp('Correct magnetic deviation')
for ii = 1 : length(mag_dev)
    [u(:,ii),v(:,ii)] = uvrot(u2(:,ii), v2(:,ii), -mag_dev(ii));   
end

%% Correct percent good: Exclude data with percent good below prct_good
disp('****')
prct_good           = input('Determine percent good threshold (generally 20): ');
pg                  = squeeze(raw.pg(:,4,first:last));  
igap                = find(pg<prct_good);           % Exclude data with percent good below prct_good
u(igap)             = NaN;
v(igap)             = NaN;
w(igap)             = NaN;
vel_err(igap)       = NaN;

%% Calculate depth of each bin:
dpt                 = sw_dpth(press,mooring.lat)';  % convert pressure to depth, press needs to have dimension (n x 1)
dpt1                = repmat(dpt,nbin,1);
binmat              = repmat((1:nbin)',1,length(dpt1));

% If ADCP is upward-looking a depth correction can be inferred from the surface reflection, which is done in adcp_surface_fit
disp('****')
if strcmp(adcp.direction,'up')  
    [z,dpt1,offset,xnull] = adcp_surface_fit(dpt,ea,sbins,blen,tlen,lag,blnk,nbin,fbind);
elseif strcmp(adcp.direction,'dn')
    z                     = dpt1+(binmat-0.5)*blen+blnk;
else
    error('Bin depth calculation: unknown direction!');
end

saveas(figure(1),[fpath_output,mooring.name,'_',num2str(adcp.sn),'_instr_',num2str(instr),'_','Histdiff_depth'],'fig')
saveas(figure(2),[fpath_output,mooring.name,'_',num2str(adcp.sn),'_instr_',num2str(instr),'_','Offset_depth'],'fig')
saveas(figure(3),[fpath_output,mooring.name,'_',num2str(adcp.sn),'_instr_',num2str(instr),'_','Amplitude_bins_2'],'fig')

%% Remove bad data if ADCP is looking upward
u1=u; v1=v; w1=w; vel_err1=vel_err; ea1=ea;

if strcmp(adcp.direction,'up')
    for i = 1:length(time)
        sz_dpt(i) = adcp_shadowzone(dpt1(i),raw.config.sysconfig.angle); % depending on the instrument depth and the beam angle the shadow zone, i.e. the depth below the surface which is contaminated by the surface reflection is determined
        iz(i)     = find(z(:,i)>sz_dpt(i),1,'last');
        sbin(i)   = bin(iz(i));
        
        %sbin(i)=30; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % here a manual criterion should be hard-coded if
        % adcp_check_surface (below) shows bad velocities close to the
        % surface
        u1(sbin(i)+1:end,i)       = NaN;
        v1(sbin(i)+1:end,i)       = NaN;
        w1(sbin(i)+1:end,i)       = NaN;
        vel_err1(sbin(i)+1:end,i) = NaN;
        ea1(sbin(i)+1:end,i)      = NaN;
    end

    if 1
        bins = nmedian(sbin)-4:nmedian(sbin)+2;
        adcp_check_surface(bins,u,u1,v,v1,time,bin,z); 
        % here the closest bins below the surface are plotted that are supposed to have good velocities, if there are still bad velocities a manual criterion needs to be found
    end
end
saveas(figure(4),[fpath_output,mooring.name,'_',num2str(adcp.sn),'_instr_',num2str(instr),'_','Meridional_zonal_velocity'],'fig')

%% SAVE DATA
% More meta information
%adcp.comment='';
adcp.config         = raw.config;
adcp.z_offset       = offset;
adcp.ang            = ang;
adcp.mag_dev        = rot;

% Data structure
data.u              = u1;
data.v              = v1;
data.w              = w1;
data.e              = vel_err1;
data.ea             = ea1;
data.pg             = pg;
data.time           = time;
data.z_bins         = z;
data.depth          = dpt;
data.temp           = temp;
data.sspd           = soundspeed;
data.lat            = mooring.lat;
data.lon            = mooring.lon;

disp('****')
reply_ts            = input('Do you want to calculate Target Strength? (CTD profile needed) 1/0 [0]:');
if isempty(reply_ts)
    reply_ts = 0;
end
if reply_ts == 1
    % EA0: noisefloor from ADCP electronic noise; a scalar
    disp('****')
    disp('Calculate Target Strength')
    disp('EA0=18 (noisefloor from ADCP electronic noise) ; pH=8.1; (seawater pH)')
    EA0 = 18;
    file_ctd        = input('Path to CTD .nc file [''*.nc'']:');
    if exist(file_ctd, 'file')
        ctd.depth = ncread(file_ctd, 'PRES');
        ctd.temp  = ncread(file_ctd, 'TEMP');
        ctd.sal   = ncread(file_ctd, 'PSAL');
        ctd.tempR = ones(size(z));
        ctd.salR  = ones(size(z));
        % for each ADCP bin, find salinity and temperature values
        for i_bin = 1:nbin
            [~,ind_depth]      = min(abs(z(i_bin,:)-ctd.depth(:)));
            ctd.salR(i_bin,:)  = ctd.sal(ind_depth);
            ctd.tempR(i_bin,:) = ctd.temp(ind_depth);
        end
        % Compute TS
        [sspd,~]  = meshgrid(soundspeed,1:nbin);
        data.ts   = target_strength(data.ea,EA0,ctd.salR,ctd.tempR,sspd,temp',soundspeed',adcp.config.cell,adcp.config.blank,adcp.config.sysconfig.angle,freq);
        figure;
        colormap jet;
        pcolor(data.ts.ts);
        shading flat;
        title('Target Strength'); colorbar;
        ylabel('Bins');
        xlabel('Time index');
        saveas(gcf,[fpath_output,mooring.name,'_',num2str(adcp.sn),'_instr_',num2str(instr),'_','target_strength'],'fig')
    else
        disp('CTD file doesn''t exist')
        reply_ts = 0;
    end
end
       
%% Save raw data
disp('****')
disp('Saving raw data')
save([fpath_output, mooring.name '_' num2str(adcp.sn) '_instr_' sprintf('%02d',instr) '.mat'],'adcp','mooring','data','raw','-v7.3');

%% Interpolate data on a regular vertical grid
Z                   = fliplr(blen/2:blen:max(z(:))+blen);
Zmax                = max(Z);
u_interp            = NaN(length(time),length(Z));
v_interp            = NaN(length(time),length(Z));
if reply_ts == 1
    ts_interp       = NaN(length(time),length(Z));
end

for i=1:length(time)
    % indice correspondant sur la grille finale Z
    ind = round((Zmax-z(1,i))/blen)+1;
    % filling the grid
    npts = min([length(Z)+ind+1 length(bin)]);
    u_interp(i,ind:ind+npts-1) = u1(1:npts,i);
    v_interp(i,ind:ind+npts-1) = v1(1:npts,i);
    if reply_ts == 1
        ts_interp(i,ind:ind+npts-1) = data.ts.ts(1:npts,i);
    end
end

if reply_ts == 0
    ts_interp       = NaN;
end
%% Horizontal interpolation, filtering and subsampling
disp('****')
disp('Filtering using 40h Hamming filter (other disponible filter: FFT filter, tide model)'),
[uintfilt,vintfilt,tsintfilt,inttim,utid_baro,vtid_baro] = adcp_filt_sub(data,u_interp',v_interp',ts_interp',1:length(Z),reply_ts);
saveas(figure(5),[fpath_output,mooring.name,'_',num2str(adcp.sn),'_instr_',num2str(instr),'_','data_raw_filt_subsampled_1'],'fig')
saveas(figure(6),[fpath_output,mooring.name,'_',num2str(adcp.sn),'_instr_',num2str(instr),'_','data_raw_filt_subsampled_2'],'fig')

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
save([fpath_output, mooring.name '_' num2str(adcp.sn) '_instr_' sprintf('%02d',instr) '_int_filt_sub.mat'],'adcp','mooring','data');

%% Figure
niv_u               = (-1:0.05:1);
niv_v               = (-1:0.05:1);

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
ylim([0,adcp.instr_depth]);
%change figure label in HH:MM
gregtick;
title({[mooring.name, ' - ZONAL VELOCITY - RDI ',num2str(freq),' kHz']});

%v
subplot(2,1,2);
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

graph_name = [fpath_output, mooring.name '_U_V_int_filt_sub'];
set(hf,'Units','Inches');
pos        = get(hf,'Position');
set(hf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(hf,graph_name,'-dpdf','-r300');

if reply_ts == 1
hf=figure('position', [0, 0, 1400, 500]);
colormap jet
[C,h] = contourf(inttim,Z(bin_start:bin_end),tsintfilt(bin_start:bin_end,:)); set(h,'LineColor','none');
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
print(hf,graph_name,'-dpdf','-r300');
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