function adcp_check_surface(bins,u,u1,v,v1,time,bin,z)

%c=colormap(jet(length(bins)));

figure(4);
orient tall;
subplot 321;
imagesc(time,bin(bins),u(bins,:),[-1 1]);
set(gca,'YTick',bins);
hold on;
gregtick;
title('Meridional velocity before correction');

subplot 323;
imagesc(time,bin(bins),u1(bins,:),[-1 1]);
set(gca,'YTick',bins);
hold on;
gregtick;
title('Meridional velocity after correction');

subplot 322;
imagesc(time,bin(bins),v(bins,:),[-1 1]);
set(gca,'YTick',bins);
hold on;
gregtick;
title('Zonal velocity after correction');

subplot 324;
imagesc(time,bin(bins),v1(bins,:),[-1 1]);
set(gca,'YTick',bins);
hold on;
gregtick;
title('Zonal velocity after correction');

subplot 325;
for i=1:length(bins);
    plot(time,hammfilter_nodec(z(bins(i),:),73),'color','k');
    hold on;
end
set(gca,'XLim',[time(1) time(end)]);
gregtick;
title('Bin depth after correction');

subplot 326;
for i=1:length(bins);
    plot(time,hammfilter_nodec(z(bins(i),:),73),'color','k');
    hold on;
end
set(gca,'XLim',[time(1) time(end)]);
gregtick;
title('Bin depth after correction');
