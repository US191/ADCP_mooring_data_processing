function out = target_strength(EA,EA0,S,T,sspd,tempx,sspdx,binlength,blank,beamangle,freq)
% function out = targetstrength(EA,EA0,S,T,sspd,tempx,sspdx,binlength,blank,beamangle)
%
% freq : khz
% EA:  echo amplitudes in counts;         [d x t]
% EA0: noisefloor from ADCP electronic noise; a scalar
% S:   salinity on ADCP-grid;             [d x t]
% T:   in-situ temperature on ADCP-grid;  [d x t]
% sspd: soundspeed on ADCP-grid;          [d x t]
% tempx: transducer temerature e.g. from Thermosal or ADCP binary data; a vector [1 x t]
% sspdx: transducer soundspeed e.g. from Thermosal or ADCP binary data; a vector [1 x t]
% binlength: ADCP configurated length of bin; a scalar
% blank:     ADCP configurated blank length; a scalar
% beamangle: ADCP beam angle to vertical; a scalar
%
% units: EA, EA0:     counts (rawdata-output)
%        T, tempx:    degrees Celsius
%        sspd, sspdx: m/s
%        beamangle:   degrees
%        bin, blank:  m
%
% out is a structure with
% out.ts: target strength in dB (relative)
% out.csspd: 'cumulative' soundspeed = effective sspd for this bin
% out.R: slant range along beam
% out.alpha: absorption coefficient = average of water column until this bin
%
% calculations following Plimpton, Freitag, McPhaden
% 'Processing of subsurface ADCP data in the equatorial Pacific'
%
% last changed 11/02/18 TF

pH          = 8.1;  % seawater pH

szm         = size(EA);
nbins       = szm(1);

dist        = ([1:nbins]'-0.5)*binlength+blank;
si          = [dist(1);diff(dist)];
tges        = cumsum(repmat(si,1,szm(2))./sspd);
csspd       = repmat(dist,1,szm(2))./tges;
out.csspd   = csspd;

Rconst      = (blank+0.31/2+([1:nbins]'-0.5)*binlength+binlength/4)/cos(beamangle*pi/180);
R           = repmat(Rconst,1,szm(2)).*csspd./repmat(sspdx,nbins,1);
out.R       = R;

%% Calculation of sound absorption (Appendix 2 Plimpton et al 2004)

% Pure water contribution:
A3          = 4.937e-4-2.59e-5*T+9.11e-7*T.^2-1.5e-8*T.^3;
A31         = 3.964e-4-1.146e-5*T+1.45e-7*T.^2-6.5e-10*T.^3;
A3(T>20)    = A31(T>20);
zz          = R*cos(beamangle*pi/180);
P3          = 1-3.83e-5*zz+4.9e-10*zz.^2;

% MgSO4 contribution:
A2          = 21.44*S.*(1+0.025*T)./sspd;
P2          = 1-1.37e-4*zz+6.2e-9*zz.^2;
f2          = 8.17*(10.^(8-1990./(273.16+T)))./(1+0.0018*(S-35));

% Boric acid contribution:
A1          = 8.86*10^(0.78*pH-5)./sspd;
f1          = 2.8*sqrt(S/35).*10.^(4-1245./(273.16+T));

alpha       = (A3.*P3*freq^2+A2.*P2.*f2.*freq^2./(freq^2+f2.^2)+A1.*f1.*freq^2./(freq^2+f1.^2))/1000; % dB/m
alpham      = cumsum(alpha)./repmat([1:nbins]',1,szm(2));
out.alpha   = alpham;

ampcor      = EA-EA0;
ampcor(ampcor<0) = 0;
out.ts      = 10*log10((repmat(tempx,nbins,1)+273.16)./sspd)+10*log10(10.^(127.3./(repmat(tempx,nbins,1)+273.16).*ampcor/10)-1)+2*alpham.*R+20*log10(R);
