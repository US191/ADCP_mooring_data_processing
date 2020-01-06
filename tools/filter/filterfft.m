%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot filtering data                                                     %
% CNRS / Météo-France / Coriolis  --- AtlantOS                            %
% Date: 14/09/2016  ---  Auteur: Pierre Rousselot                         %
% fft_filter: Patrick Martineau                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tempf(end+1:length(tempf)*2,:)=tempf;
% positions.timef(end+1:length(positions.timef)*2,:)=positions.timef;
% depth(end+1:length(depth)*2,:)=depth;
% positions.timef = linspace(1,15940,15940);
% for i = 2:17
%     positions.timef(i,:)=positions.timef(1,:);
% end

[reconstructed2]=fft_filter(tempf,1,[0 24]);

figure
plot(positions.timef,reconstructed)
datetick('x','mm/yy','keepticks')
xlabel('Time')
ylabel('Temperature [degC]')
%axis([min(positions.timef(:,1)) max(positions.timef(:,1)) 11 22])
grid on
title('Filtering process diff (<1week) - (<1day)')

% figure;
% contourf(positions.timef', -depth, reconstructed)

figure
subplot(4,1,1)
contourf(vec_mean,-depth_mean,temp_mean)
axis([7.34313e5 7.34318e5 -83 0])
datetick('x','dd/mm','keeplimits')
xlabel('Time')
ylabel('Depth [m]')
cc=colorbar;
ylabel(cc,'Temperature [{\circ}C]')
title('Daily mean')

subplot(4,1,2)
contourf(positions.timef,-depth,tempf)
axis([7.34313e5 7.34318e5 -83 0])
datetick('x','dd/mm','keeplimits')
xlabel('Time')
ylabel('Depth [m]')
cc=colorbar;
ylabel(cc,'Temperature [{\circ}C]')
title('Raw Data')

subplot(4,1,3)
contourf(positions.timef,-depth,reconstructed)
axis([7.34313e5 7.34318e5 -83 0])
datetick('x','dd/mm','keeplimits')
xlabel('Time')
ylabel('Depth [m]')
cc=colorbar;
ylabel(cc,'Temperature [{\circ}C]')
title('Filtered Data (>24h)')

subplot(4,1,4)
contourf(positions.timef,-depth,reconstructed2)
axis([7.34313e5 7.34318e5 -83 0])
xlabel('Time')
ylabel('Depth [m]')
cc=colorbar;
ylabel(cc,'Temperature Anomaly [{\circ}C]')
title('Diurnal Temperature Anomaly (<24h)')



