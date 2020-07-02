
wind = dlmread('/home/proussel/Documents/OUTILS/ADCP/ADCP_mooring_data_processing/FR28/w0n0e_dy.ascii');
x = (1:5:684)';
y = zeros(size(wind(1:5:end,3)));
u = wind(1:5:end,3);
v = wind(1:5:end,4);
time = strcat(num2str(wind(1:5:end,1)),num2str(wind(1:5:end,2)));
for ii = 1 :length(time)
    ti(ii) = julian(str2double(time(ii,1:4)),str2double(time(ii,5:6)),str2double((time(ii,7:8))));
end
%time=datenum(time,'yyyymmddHHMM');
quiver(ti',y,u,v)
%datetick('x')
gregtick