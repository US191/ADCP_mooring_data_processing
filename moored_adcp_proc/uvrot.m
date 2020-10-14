function [ru,rv]=uvrot(u,v,ang)
%function [ru,rv]=uvrot(u,v,ang)

ang=ang*pi/180;
ru=u.*cos(ang)-v.*sin(ang);
rv=u.*sin(ang)+v.*cos(ang);


























