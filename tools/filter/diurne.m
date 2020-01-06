%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Biais capteurs - zoom cycle diurne                                      %
% CNRS / Météo-France / Coriolis  --- AtlantOS                            %
% Date: 01/08/2016  ---  Auteur: Pierre Rousselot                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%addpath('
date = datestr(positions.timef(:,1));

% Profil fin de nuit
time =datestr(positions.timefinterp(:,1),'HH:MM');
j=1;
jj=1;
for ii = 2:length(time(:,1))
    %end_night(ii)=strcmp(time(ii,:),'03:00');
    if strcmp(time(ii,:),'03:00')
        end_night(j)=ii;
        j=j+1;
        if strcmp(time(ii+3,:),'06:00')
            end_night_bis(jj)=ii+3;
            jj=jj+1;
        else
            if strcmp(time(ii+2,:),'05:00')
                end_night_bis(jj)=ii+2;
                jj=jj+1;              
            else
                if strcmp(time(ii+1,:),'04:00')
                    end_night_bis(jj)=ii+1;
                    jj=jj+1;    
                else
                    if strcmp(time(ii,:),'03:00')
                        end_night_bis(jj)=ii;
                        jj=jj+1; 
                    end
                end
            end
        end
    end

end

for ii= 1:331
    mean_temp(ii,:) = mean(tempf(end_night(ii):end_night_bis(ii),:),1);
    mean_depth(ii,:) = mean(depth(end_night(ii):end_night_bis(ii),:),1);
    mean_time(ii,:) = mean(positions.timefinterp(end_night(ii):end_night_bis(ii),:),1);
end

% %Conv filtering
% for ii=1:length(tempf(1,:))
%     ttt(:,ii)=conv(tempf(end_night,ii),ones(1,3)/3,'valid');
% end
% 
% %Median filtering
% %tt=medfilt1(tempf(end_night,:));
% hh=hann(24);
% for ii=1:length(tempf(1,:))
%     tt(:,ii)=conv(tempf(end_night,ii),hh,'valid');
% end
% 
% 
% %Filter
% %tttt=filter(ones(1,3)/3,1,(tempf(end_night,:)));
% 
% %Zero-Phase filtering
% ttttt=filtfilt(ones(1,3)/3,1,(tempf(end_night,:)));
% 
% %Savitky-Golay fitering
% tttttt=sgolayfilt(tempf(:,:),3,7);
% 
% %LPO-ecart median
% % for ii=1:length(tempf(1,:))
% %     [par(:,ii),info_rejet_med]=clean_median(tempf(:,ii),6,2.8,[0.05 1],1,NaN);
% % end
% 
% figure(1);
% title('End-Night - Different filter')
% subplot(5,1,1)
% plot(tempf(end_night,:))
% subplot(5,1,2)
% plot(tt)
% subplot(5,1,3)
% plot(ttt)
% subplot(5,1,4)
% plot(ttttt)
% subplot(5,1,5)
% plot(tttttt)

figure
contourf(mean_time, -mean_depth, mean_temp);



figure(2);
datestr(positions.timefinterp(end_night,1));
contourf(positions.timefinterp(end_night,:), -depth_ma(end_night,:), yyyy(:,end_night)')
datetick('x','mm/yy','keepticks')
cc = colorbar;
xlabel('Time')
ylabel('Depth [m]')
ylabel(cc,'Temperature [{\circ}C]');
caxis([min(min(tempf_interp)) max(max(tempf_interp))])
grid on
title('End-Night Mean (3h-6h) Time Series')

% TT = 18;
% %Zoom cycle diurne - biais capteur (fin de nuit)
% figure(3);
% plot(positions.timef(TT:TT+700,:),tempf(TT:TT+700,:))
% datetick('x','HH dd/mm','keepticks')
% grid on
% j1m = mean(tempf(46:49,:));
% j2m = mean(tempf(69:72,:));
% diff_mean = j1m-j2m
% title('Zoom 2 days')


%Biais intercapteur
% figure(4);
% 
% if options.marisonde
%     %plot(positions.timef(500+24*3:500+24*6,:),tempf(500+24*3:500+24*6,:))
%     plot(positions.timef(TT+24*3:TT+24*6,:),tempf(TT+24*3:TT+24*6,:))
% else
%     plot(positions.timef(5100+24*3:5100+24*6,:),tempf(5100+24*3:5100+24*6,:))
% end
% datetick('x','HH dd/mm','keepticks')
% grid on
% title('Zoom Homogeneous Period')


% figure(5);
% if options.marisonde
%     %plot(tempf(500,:))
%     plot(tempf(TT,:))
%     hold on
%     %plot(mean(tempf(500,:))*ones(1,17),'r--')
%     plot(mean(tempf(TT,:))*ones(1,17),'r--')
%     %max_ecart_intersensor = range(tempf(500,:))  
%     max_ecart_intersensor = range(tempf(TT,:)) 
% else
%     plot(tempf(5178+24,:))
%     hold on
%     plot(mean(tempf(5178+24,:))*ones(1,17),'r--')
%     max_ecart_intersensor = range(tempf(5178,:))
% end
% hold off
% grid on
% view([90 90])
% title('Inter-sensor Bias')
% 
% %[maxbias_intersensor,ind_maxbias_intersensor]=max(std(tempf(4934:5280,:),[],2))
% [maxbias_intersensor,ind_maxbias_intersensor]=max(std(tempf(TT:100,:),[],2))
% 
% % mm = mean(tempf(5100+24*3:5100+24*6,:));
% % stdmm = std(mm);
% mm = mean(tempf(TT+24*3:TT+24*6,:));
% stdmm = std(mm);
% 
% %Derive dans le temps??
% j=1;
% %for ii=4950:7400 %1:length(time(:,1))
% for ii=TT:TT+1100 %1:length(time(:,1))
%     %end_night(ii)=strcmp(time(ii,:),'03:00');
%     if strcmp(time(ii,:),'03:00')
%         end_night2(j)=ii;
%         j=j+1;
%     end
% end
% mmean = nanmean(tempf(end_night2,:),2);
% for ii = 1:length(tempf(1,:))
% diff(:,ii) = tempf(end_night2,ii)-mmean;
% end
% mmmmm = nanmean(diff);
% for ii=1:length(end_night2(1,:))
% te(ii,:)=tempf(end_night2(1,ii),:)-mmmmm;
% end