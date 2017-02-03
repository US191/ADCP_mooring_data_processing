function z_shadowzone=adcp_shadowzone(adcp_dpt,beamangle);

% z_shadowzone=adcp_shadowzone(adcp_dpt,beamangle);
%
% Computes depth below the surface/above ground that is invalid due to surface/bottom reflections
% Input: - Depth below surface/above ground of ADCP [m]
%        - Angle of beams [deg]
%
% R. Kopte, 2014/11/20

z_shadowzone=adcp_dpt*(1-cosd(beamangle));
