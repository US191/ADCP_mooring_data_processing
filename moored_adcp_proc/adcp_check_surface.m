function adcp_check_surface(bins,u,u1,v,v1,corr,corr1,time,bin,z,z1,sz_dpt)

c=colormap(jet(length(bins)));
figure(8);
colormap jet
orient tall;
subplot 331;
imagesc(time,bin(bins),u(bins,:),[-1 1]);
set(gca,'YTick',bins);
set(gca,'YDir','normal')
ylabel('Bins');
hold on;
gregtick;
title('Zonal velocity before correction');

subplot 334;
imagesc(time,bin(bins),u1(bins,:),[-1 1]);
set(gca,'YTick',bins);
set(gca,'YDir','normal')
hold on;
ylabel('Bins');
gregtick;
title('Zonal velocity after correction');

subplot 332;
imagesc(time,bin(bins),v(bins,:),[-1 1]);
set(gca,'YTick',bins);
set(gca,'YDir','normal')
hold on;
ylabel('Bins');
gregtick;
title('Meridional velocity before correction');

subplot 335;
imagesc(time,bin(bins),v1(bins,:),[-1 1]);
set(gca,'YTick',bins);
set(gca,'YDir','normal')
hold on;
ylabel('Bins');
gregtick;
title('Meridional velocity after correction');

subplot 333;
imagesc(time,bin(bins),corr(bins,:));
set(gca,'YTick',bins);
set(gca,'YDir','normal')
hold on;
ylabel('Bins');
gregtick;
title('Correlation before correction');

subplot 336;
imagesc(time,bin(bins),corr1(bins,:));
set(gca,'YTick',bins);
set(gca,'YDir','normal')
hold on;
ylabel('Bins');
gregtick;
title('Correlation after correction');

subplot 337;
for i=1:length(bins)
    plot(time,-hammfilter_nodec(z(bins(i),:),73),'color','k');
    hold on;
    plot(time,-hammfilter_nodec(z1(bins(i),:),73),'color','g');
end
plot(time,-sz_dpt,'--r')
set(gca,'XLim',[time(1) time(end)]);
gregtick;
ylabel('Depth');
title('Bin depth after correction');

subplot 338;
for i=1:length(bins)
    plot(time,-hammfilter_nodec(z(bins(i),:),73),'color','k');
    hold on;
    plot(time,-hammfilter_nodec(z1(bins(i),:),73),'color','g');
end
plot(time,-sz_dpt,'--r')
set(gca,'XLim',[time(1) time(end)]);
gregtick;
ylabel('Depth');
title('Bin depth after correction');

subplot 339;
for i=1:length(bins)
    plot(time,-hammfilter_nodec(z(bins(i),:),73),'color','k');
    hold on;
    plot(time,-hammfilter_nodec(z1(bins(i),:),73),'color','g');
end
plot(time,-sz_dpt,'--r')
set(gca,'XLim',[time(1) time(end)]);
gregtick;
ylabel('Depth');
title('Bin depth after correction');
