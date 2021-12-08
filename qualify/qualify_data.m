%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison ADCP Data Mooring with SADCP/LADCP/DVL
% Autor: Pierre Rousselot
% Date:  03/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; 
addpath(genpath('moored_adcp_proc'))
addpath(genpath('tools'))

%% Load Mooring data
Moorfile  = 'FR31/0-0/0W0N_22545_instr_01.mat';
Moorfile2 = 'FR31/0-0/0W0N_22545_instr_01_int_filt_sub.mat';
% Moorfile  = 'FR31/10W/10W0N_24629_instr_01.mat';
% Moorfile2 = 'FR31/10W/10W0N_24629_instr_01_int_filt_sub.mat';

%FR31
OS38file  = '/media/irdcampagnes/PIRATA/PIRATA-FR31/data-final/SADCP/OS38/ncc/FR31-OS38_osite_final_fhv1.nc';
OS150file = '/media/irdcampagnes/PIRATA/PIRATA-FR31/data-final/SADCP/OS150/ncc/FR31-OS150_osite_final_fhv1.nc';
DVLfile   = '/media/irdcampagnes/PIRATA/PIRATA-FR31/data-final/SADCP/DVL600/ncc/FR31-DVL600_osite_final_fhv1.nc';
LADCPfil1 = '/media/irdcampagnes/PIRATA/PIRATA-FR31/data-adjusted/LADCP/profiles/FR31_006.mat';
LADCPfil2 = '/media/irdcampagnes/PIRATA/PIRATA-FR31/data-adjusted/LADCP/profiles/FR31_031.mat';
% %
% OS38file  = '/media/irdcampagnes/PIRATA/PIRATA-FR30/data-final/SADCP/OS38/ncc/FR30-OS38_osite_mat20_fhv1_corr_final.nc';
% OS150file = '/media/irdcampagnes/PIRATA/PIRATA-FR30/data-final/SADCP/OS150/ncc/FR30-OS150_osite_mat20_fhv1_corr_final.nc';
% DVLfile   = '/media/irdcampagnes/PIRATA/PIRATA-FR30/data-final/SADCP/LOCH/ncc/FR30-DVL600_osite_mat20_fhv1_corr_final.nc';
% LADCPfil1 = '/media/irdcampagnes/PIRATA/PIRATA-FR30/data-adjusted/LADCP/profiles/FR30_009.mat';
% LADCPfil2 = '/media/irdcampagnes/PIRATA/PIRATA-FR30/data-adjusted/LADCP/profiles/FR30_010.mat';
% % 
% OS38file  = '/media/irdcampagnes/PIRATA/PIRATA-FR29/data-final/SADCP/OS38/data/LTA/PIRATA-FR29-OS38_LTA_osite_fhv1_final.nc';
% OS150file = '/media/irdcampagnes/PIRATA/PIRATA-FR29/data-final/SADCP/OS150/data/LTA/PIRATA-FR29-OS150_LTA_osite_fhv1_final.nc';
% DVLfile   = '/media/irdcampagnes/PIRATA/PIRATA-FR29/data-final/SADCP/DVL600/data/LTA/PIRATA-FR29-DVL600_LTA_osite_fhv1_final.nc';
% LADCPfil1 = '/media/irdcampagnes/PIRATA/PIRATA-FR29/data-adjusted/LADCP/profiles/FR29_00003.mat';
% LADCPfil2 = '/media/irdcampagnes/PIRATA/PIRATA-FR29/data-adjusted/LADCP/profiles/FR29_00045.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load dataset
%% Load mooring Data
load(Moorfile)
if contains(Moorfile, 'int_filt_sub')
    u_mooring     = data.uintfilt;
    v_mooring     = data.vintfilt;
    t_mooring     = data.inttim;
    z_mooring     = data.Z;
elseif contains(Moorfile, 'int_filt')
    u_mooring     = uifilt;
    v_mooring     = vifilt;
    t_mooring     = data.time;
    z_mooring     = data.Z;
else
    u_mooring     = data.u;
    v_mooring     = data.v;  
    t_mooring     = data.time;
    z_mooring     = data.z_bins;
end

load(Moorfile2)
u_mooring2     = data.uintfilt;
v_mooring2     = data.vintfilt;
t_mooring2     = data.inttim;
z_mooring2     = data.Z;
    
%% Load DVL Data
dvl.u    = ncread(DVLfile, 'UVEL_ADCP_CORTIDE');
dvl.v    = ncread(DVLfile, 'VVEL_ADCP_CORTIDE');
dvl.t    = ncread(DVLfile, 'JULD');
dvl.d    = ncread(DVLfile, 'DEPH');
dvl.lat  = ncread(DVLfile, 'LATITUDE'); 
dvl.lon  = ncread(DVLfile, 'LONGITUDE');

%% Load OS150 Data
os.u     = ncread(OS150file, 'UVEL_ADCP_CORTIDE');
os.v     = ncread(OS150file, 'VVEL_ADCP_CORTIDE');
os.t     = ncread(OS150file, 'JULD');
os.d     = ncread(OS150file, 'DEPH');
os.lat   = ncread(OS150file, 'LATITUDE'); 
os.lon   = ncread(OS150file, 'LONGITUDE');

%% Load OS38 Data
os38.u   = ncread(OS38file, 'UVEL_ADCP_CORTIDE');
os38.v   = ncread(OS38file, 'VVEL_ADCP_CORTIDE');
os38.t   = ncread(OS38file, 'JULD');
os38.d   = ncread(OS38file, 'DEPH');
os38.lat = ncread(OS38file, 'LATITUDE'); 
os38.lon = ncread(OS38file, 'LONGITUDE');

%% Load LADCP Data
load(LADCPfil1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Research Date - reference LADCP
h_mooring     = datenum(julian(t_mooring(:)));
h_mooring2    = datenum(julian(t_mooring2(:)));

h_ladcp       = datenum(dr.date);
[~,plt_m]     = min(abs((h_mooring)-h_ladcp));
[~,plt_m2]    = min(abs((h_mooring2)-h_ladcp));

h_os          = os.t + datenum('1950-01-01');
[~,plt_os]    = min(abs((h_os)-h_ladcp));

h_os38        = os38.t + datenum('1950-01-01');
[~,plt_os38]  = min(abs((h_os38)-h_ladcp));

h_dvl         = dvl.t + datenum('1950-01-01');
[~,plt_dvl]   = min(abs((h_dvl)-h_ladcp));

if contains(Moorfile, 'int_filt')
    z_mooring_ext = z_mooring;
else
    z_mooring_ext = z_mooring(:, plt_m);
end
z_mooring_ext2 = z_mooring2;

%% Calc distance
d_station     = dist(data.lat,data.lon,dr.lat,dr.lon);
d_station     = d_station * 0.539957;

%% Plot Before
% plt_dvl = plt_dvl + 1;
% plt_os  = plt_os-1;

figure
subplot(1,2,1)
plot(dr.u, -dr.z, '--g')
hold on; 
plot(os.u(:,plt_os), os.d,'r')
plot(os38.u(:,plt_os38), os38.d,'c')
plot(dvl.u(:,plt_dvl), dvl.d,'k')
plot(u_mooring(:,plt_m), -z_mooring_ext,'LineWidth',2,'Color','b')
plot(u_mooring2(:,plt_m2), -z_mooring_ext2,'--','LineWidth',1.5,'Color','b')

hold off
xlabel('U [m.s^{-1}]')
ylabel('Depth [m]')
grid on
ylim([-300 0])
date_dvl   = (h_dvl(plt_dvl)-h_mooring(plt_m))*(24*60);
date_os    = (h_os(plt_os)-h_mooring(plt_m))*(24*60);
date_os38  = (h_os38(plt_os38)-h_mooring(plt_m))*(24*60);
date_ladcp = (h_ladcp-h_mooring(plt_m))*(24*60);
if date_os<0
        legend(['LADCP / ' num2str(round(date_ladcp)) 'min'], ['OS150 / ' num2str(round(date_os)) 'min'], ['OS38 / ' num2str(round(date_os38)) 'min'], ['DVL / ' num2str(round(date_dvl)) 'min'], ['ADCP Mooring / ' datestr(h_mooring(plt_m))], 'ADCP filtered','location', 'southeast')
else
        legend(['LADCP / +' num2str(round(date_ladcp)) 'min'], ['OS150 / +' num2str(round(date_os)) 'min'], ['OS38 / +' num2str(round(date_os38)) 'min'], ['DVL / +' num2str(round(date_dvl)) 'min'], ['ADCP Mooring / ' datestr(h_mooring(plt_m))],'ADCP filtered', 'location', 'southeast')
end
title(['ADCP Mooring measurement vs other current measurement sources [distance = ' num2str(round(d_station*10)/10) 'NM]'])

subplot(1,2,2)
plot(dr.v, -dr.z, '--g')
hold on; 
plot(os.v(:,plt_os), os.d,'r')
plot(os38.v(:,plt_os38), os38.d,'c')
plot(dvl.v(:,plt_dvl), dvl.d,'k')
plot(v_mooring(:,plt_m), -z_mooring_ext,'LineWidth',2,'Color','b')
plot(v_mooring2(:,plt_m2), -z_mooring_ext2,'--','LineWidth',1.5,'Color','b')
hold off
xlabel('V [m.s^{-1}]')
ylabel('Depth [m]')
grid on
ylim([-300 0])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load LADCP Data
load(LADCPfil2)

%% Research Date
h_ladcp     = datenum(dr.date);
[~,plt_m]   = min(abs((h_mooring)-h_ladcp));
[~,plt_m2]  = min(abs((h_mooring2)-h_ladcp));

[~,plt_os]  = min(abs((h_os)-h_ladcp));

[~,plt_os38]  = min(abs((h_os38)-h_ladcp));

[~,plt_dvl] = min(abs((h_dvl)-h_ladcp));

if contains(Moorfile, 'int_filt')
    z_mooring_ext = z_mooring;
else
    z_mooring_ext = z_mooring(:, plt_m);
end

%% Calc distance
d_station     = dist(data.lat,data.lon,dr.lat,dr.lon);
d_station     = d_station * 0.539957;

%% Plot After
figure
subplot(1,2,1)
plot(dr.u, -dr.z, '--g')
hold on; 
plot(os.u(:,plt_os), os.d,'r')
plot(os38.u(:,plt_os38), os38.d,'c')
plot(dvl.u(:,plt_dvl), dvl.d,'k')
plot(u_mooring(:,plt_m), -z_mooring_ext,'LineWidth',2,'Color','b')   %data.Z(:,plt_m)
plot(u_mooring2(:,plt_m2), -z_mooring_ext2,'--','LineWidth',1.5,'Color','b')
hold off
xlabel('U [m.s^{-1}]')
ylabel('Depth [m]')
grid on
ylim([-300 0])
legend(['LADCP / ' datestr(h_ladcp)], ['OS150 / ' datestr(h_os(plt_os))], ['OS38 / ' datestr(h_os38(plt_os38))], ['DVL / ' datestr(h_dvl(plt_dvl))], ['ADCP Mooring / ' datestr(h_mooring(plt_m))],'ADCP filtered','location', 'southeast')
date_dvl   = (h_dvl(plt_dvl)-h_mooring(plt_m))*(24*60);
date_os    = (h_os(plt_os)-h_mooring(plt_m))*(24*60);
date_os38  = (h_os38(plt_os38)-h_mooring(plt_m))*(24*60);
date_ladcp = (h_ladcp-h_mooring(plt_m))*(24*60);
if date_os<0
        legend(['LADCP / ' num2str(round(date_ladcp)) 'min'], ['OS150 / ' num2str(round(date_os)) 'min'], ['OS38 / ' num2str(round(date_os38)) 'min'], ['DVL / ' num2str(round(date_dvl)) 'min'], ['ADCP Mooring / ' datestr(h_mooring(plt_m))], 'ADCP filtered','location', 'southeast')
else
        legend(['LADCP / +' num2str(round(date_ladcp)) 'min'], ['OS150 / +' num2str(round(date_os)) 'min'], ['OS38 / +' num2str(round(date_os38)) 'min'], ['DVL / +' num2str(round(date_dvl)) 'min'], ['ADCP Mooring / ' datestr(h_mooring(plt_m))],'ADCP filtered','location', 'southeast')
end
title(['ADCP Mooring measurement vs other current measurement sources [distance = ' num2str(round(d_station*10)/10) 'NM]'])

subplot(1,2,2)
plot(dr.v, -dr.z, '--g')
hold on; 
plot(os.v(:,plt_os), os.d,'r')
plot(os38.v(:,plt_os38), os38.d,'c')
plot(dvl.v(:,plt_dvl), dvl.d,'k')
plot(v_mooring(:,plt_m), -z_mooring_ext,'LineWidth',2,'Color','b')
plot(v_mooring2(:,plt_m2), -z_mooring_ext2,'--','LineWidth',1.5,'Color','b')
hold off
xlabel('V [m.s^{-1}]')
ylabel('Depth [m]')
grid on
ylim([-300 0])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Research nearest position
d_os38          = dist(data.lat,data.lon,os38.lat,os38.lon);
[d_mo,plt_os38] = min(d_os38);
d_mo            = d_mo * 0.539957;

d_os150         = dist(data.lat,data.lon,os.lat,os.lon);
[~,plt_os]      = min(d_os150);

d_dvl           = dist(data.lat,data.lon,dvl.lat,dvl.lon);
[~,plt_dvl]     = min(d_dvl);

[~,plt_m]       = min(abs(h_os38(plt_os38)-h_mooring));
[~,plt_m2]       = min(abs(h_os38(plt_os38)-h_mooring2));
%plt_m = plt_m-3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% En dur
% plt_os38 = plt_os38 -20;
% plt_os = plt_os -20;
% plt_dvl = plt_dvl -20;
% d_mo          = dist(lat_moor,lon_moor,os38.lat(plt_os38+1),os38.lon(plt_os38+1));
% d_mo            = d_mo * 0.539957;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if contains(Moorfile, 'int_filt')
    z_mooring_ext = z_mooring;
else
    z_mooring_ext = z_mooring(:, plt_m);
end

%% Plot 
figure
subplot(1,2,1)
hold on; 
plot(os.u(:,plt_os), os.d,'r')
plot(os38.u(:,plt_os38), os38.d,'c')
plot(dvl.u(:,plt_dvl), dvl.d,'k')
plot(u_mooring(:,plt_m), -z_mooring_ext,'LineWidth',2,'Color','b')
plot(u_mooring2(:,plt_m2), -z_mooring_ext2,'--','LineWidth',1.5,'Color','b')
hold off
xlabel('U [m.s^{-1}]')
ylabel('Depth [m]')
grid on
ylim([-300 0])
date_dvl  = (h_dvl(plt_dvl)-h_mooring(plt_m))*(24*60);
date_os   = (h_os(plt_os)-h_mooring(plt_m))*(24*60);
date_os38 = (h_os38(plt_os38)-h_mooring(plt_m))*(24*60);
if date_os<0
        legend(['OS150 / ' num2str(round(date_os)) 'min'], ['OS38 / ' num2str(round(date_os38)) 'min'], ['DVL / ' num2str(round(date_dvl)) 'min'], ['ADCP Mooring / ' datestr(h_mooring(plt_m))], 'ADCP filtered','location', 'southeast')  
else
        legend(['OS150 / +' num2str(round(date_os)) 'min'],['OS38 / +' num2str(round(date_os38)) 'min'], ['DVL / +' num2str(round(date_dvl)) 'min'], ['ADCP Mooring / ' datestr(h_mooring(plt_m))], 'ADCP filtered','location', 'southeast')
end
title(['ADCP Mooring measurement vs other current measurement sources [distance = ' num2str(round(d_mo*10)/10) 'NM]'])

subplot(1,2,2)
hold on; 
plot(os.v(:,plt_os), os.d,'r')
plot(os38.v(:,plt_os38), os38.d,'c')
plot(dvl.v(:,plt_dvl), dvl.d,'k')
plot(v_mooring(:,plt_m), -z_mooring_ext,'LineWidth',2,'Color','b')
plot(v_mooring2(:,plt_m2), -z_mooring_ext2,'--','LineWidth',1.5,'Color','b')
hold off
xlabel('V [m.s^{-1}]')
ylabel('Depth [m]')
grid on
ylim([-300 0])