function gregtick(sep)
% function gregtick(sep)
%
% GEOMAR SVN $Id: gregtick.m 95 2015-01-14 16:07:01Z gkrahmann@geomar.de $
%
%GREGTICK Gregorian date X-axis tick labels.
%
%  Christian Mertens, IfM Kiel
%  $Revision: 95 $ $Date: 2015-01-14 17:07:01 +0100 (Wed, 14 Jan 2015) $
%
% getuprop and setuprop exchanged for getappdata and setappdata
% last changed 11/08/15 TF

global MTICK SEP YEAR MONTH DAY HOUR MINUTE

% maximum number of tick labels
MTICK = 18;

% label seperation
if nargin < 1
  SEP = 3;
else
  SEP = sep;
end

YEAR = 365.25;
MONTH = YEAR/12;
DAY = 1;
HOUR = DAY/24;
MINUTE = HOUR/60;

units = get(gca,'Units');
set(gca,'Units','normalized')
pos = get(gca,'Position');
set(gca,'Units',units)

jd = get(gca,'XLim');
pos0 = get(0,'DefaultAxesPosition');
dtim = diff(jd)*pos0(3)/pos(3);
MTICK = ceil(MTICK/pos0(3)*pos(3));

% remove old ticks
set(gca,'XTick',[],'XTickLabel',[])
hg = getappdata(gca,'GregTick');
delete(hg)

if (dtim <= HOUR*8)
  hg = minutes(jd);
elseif (dtim > HOUR*8 & dtim <= DAY*10)
  hg = hours(jd);
elseif (dtim > DAY*10 & dtim <= MONTH*6)
  hg = days(jd);
elseif (dtim > MONTH*6 & dtim <= YEAR*6)
  hg = months(jd);
elseif (dtim > YEAR*6 & dtim <= YEAR*40)
  years(jd);
else
  decades(jd);
end

setappdata(gca,'GregTick',hg)


%-------------------------------------------------------------------------------

function hg = minutes(jd)
%MINUTES Label minutes.

global MTICK SEP HOUR MINUTE

smonth = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct', ...
          'Nov','Dec'};
fontsize = get(gca,'FontSize');

ilabel = [1,2,5,10,15,20,30];
[dummy,i] = min(abs(diff(jd)/MINUTE./ilabel - MTICK));
inc = ilabel(i)*MINUTE;

[yy,mm,dd,hh] = gregorian(jd);
xtick = ceil(jd(1)/inc)*inc:inc:jd(2);
for i = 1:length(xtick)
  h = (xtick(i) - floor(xtick(i)))/HOUR;
  m = round((h - floor(h))*60);
  if m >= 60
    m = 0;
    h = ceil(h);
  end
  xticklabel{i} = sprintf('%2.2d:%2.2d',floor(h),m);
end

set(gca,'XTick',xtick,'XTickLabel',xticklabel)

% label days
yt = yloc(SEP,1);
[yy,mm,dd,hh] = gregorian(xtick);
i = find(dd == dd(1));
hg(1) = text(median(xtick(i)),yt,sprintf('%s %d',smonth{mm(1)},dd(1)), ...
             'FontSize',fontsize,'HorizontalAlignment','center', ...
             'VerticalAlignment','bottom');
if floor(jd(2)) > floor(jd(1))
  i = find(dd == dd(end));
  hg(2) = text(median(xtick(i)),yt,sprintf('%s %d',smonth{mm(end)},dd(end)), ...
               'FontSize',fontsize,'HorizontalAlignment','center', ...
               'VerticalAlignment','bottom');
end


%-------------------------------------------------------------------------------

function hg = hours(jd)
%HOURS Label hours.

global MTICK SEP DAY HOUR

smonth = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct', ...
          'Nov','Dec'};
fontsize = get(gca,'FontSize');

ilabel = [1:4,6,8,12];
[dummy,i] = min(abs(diff(jd)/HOUR./ilabel - MTICK));
inc = ilabel(i)*HOUR;

[yy,mm,dd,hh] = gregorian(jd);
xtick = ceil(jd(1)/inc)*inc:inc:jd(2);
for i = 1:length(xtick)
  h = round((xtick(i) - floor(xtick(i)))/HOUR);
  if h >= 24
    h = 0;
  end
  xticklabel{i} = sprintf('%2.2d',h);
end

set(gca,'XTick',xtick,'XTickLabel',xticklabel)

% label days
hg = [];
yt = yloc(SEP,1);
jd1 = ceil(xtick(1)/HOUR)*HOUR:HOUR:xtick(end);
[yy,mm,dd,hh] = gregorian(jd1);
l = find(round(hh) == 12);
if length(l) > 0
  xt = jd1(l);
  for i = 1:length(l)
    t{i} = sprintf('%s %d',smonth{mm(l(i))},dd(l(i)));
  end
  hg = text(xt,yt*ones(size(xt)),t,'FontSize',fontsize, ...
            'HorizontalAlignment','center','VerticalAlignment','bottom');
end

%%% bars
%%i = find(hh == 0);
%%xt = jd1(i);
%%l = length(xt);
%%text(xt,yt*ones(1,l),repmat('|',l,1),'FontSize',fontsize, ...
%%     'HorizontalAlignment','center','VerticalAlignment','bottom')

% first day
[yy0,mm0,dd0,hh0] = gregorian(xtick(1));
if (hh0 > 12 & hh0 < 23)
  t = sprintf('%s %d',smonth{mm(1)},dd(1));
  if diff(jd) > DAY
    ht = text(jd(1),yt,t,'FontSize',fontsize, ...
              'HorizontalAlignment','left','VerticalAlignment','bottom');
  else
    i = find(dd == dd(1));
    ht = text(median(xtick(i)),yt,t,'FontSize',fontsize, ...
              'HorizontalAlignment','left','VerticalAlignment','bottom');
  end
  hg = [hg;ht];
end

% last day
[yy0,mm0,dd0,hh0] = gregorian(xtick(end));
if (hh0 < 12 & hh0 > 1)
  t = sprintf('%s %d',smonth{mm(end)},dd(end));
  if diff(jd) > DAY
    ht = text(jd(end),yt,t,'FontSize',fontsize, ...
              'HorizontalAlignment','right','VerticalAlignment','bottom');
  else
    i = find(dd == dd(end));
    ht = text(median(xtick(i)),yt,t,'FontSize',fontsize, ...
              'HorizontalAlignment','left','VerticalAlignment','bottom');
  end
  hg = [hg;ht];
end


%-------------------------------------------------------------------------------

function hg = days(jd)
%DAYS Label days.

global MTICK SEP MONTH

lmonth = {'January','February','March','April','May','June','July', ...
          'August','September','October','November','December'};
fontsize = get(gca,'FontSize');

ltab = {[1 15],[1 10 20],[1 8 16 24],[2:2:30],[1:31]};
ntab = length(ltab);
for i = 1:ntab
  ilabel(i) = MONTH/length(ltab{i});
end
[dummy,i] = min(abs(diff(jd)./ilabel - MTICK));
id = ltab{i};
nd = length(id);

[yy,mm,dd] = gregorian(jd);
nyy = diff(yy) + 1;
yy = repmat([yy(1):yy(2)],12*nd,1);
mm = repmat([1:12],nd,nyy);
dd = repmat(id',1,12*nyy);
yy = yy(:);
mm = mm(:);
dd = dd(:);
while 1
  i = find(diff(julian(yy,mm,dd)) < 1);
  if isempty(i), break, end
  yy(i) = [];
  mm(i) = [];
  dd(i) = [];
end
xtick = julian(yy,mm,dd);
i = find(xtick >= jd(1) & xtick <= jd(2));
xtick = xtick(i);
yy = yy(i);
mm = mm(i);
dd = dd(i);

for i = 1:length(dd)
  xticklabel{i} = sprintf('%d',dd(i));
end

set(gca,'XTick',xtick,'XTickLabel',xticklabel)
hg = [];

% label months
yt = yloc(SEP,1);
jd1 = xtick(1):xtick(end);
[yy,mm,dd] = gregorian(jd1);
l = find(dd == 15 & jd1 >= jd(1) & jd1 <= jd(2));
if length(l) > 0
  for i = 1:length(l)
    t{i} = sprintf('%s',lmonth{mm(l(i))});
  end
  xt = jd1(l);
  hg = text(xt,yt*ones(size(xt)),t,'FontSize',fontsize, ...
            'HorizontalAlignment','center','VerticalAlignment','bottom');
end

% first month
[yy0,mm0,dd0] = gregorian(xtick(1));
if (dd0 > 15 & dd0 < 30)
  if diff(jd) > MONTH
    ht = text(jd(1),yt,lmonth{mm(1)},'FontSize',fontsize, ...
              'HorizontalAlignment','left','VerticalAlignment','bottom');
  else
    i = find(mm == mm(1));
    ht = text(median(jd1(i)),yt,lmonth{mm(1)},'FontSize',fontsize, ...
              'HorizontalAlignment','center','VerticalAlignment','bottom');
  end
  hg = [hg;ht];
end

% last month
[yy0,mm0,dd0] = gregorian(xtick(end));
if (dd0 < 15 & dd0 > 2)
  if diff(jd) > MONTH
    ht = text(jd(2),yt,lmonth{mm(end)},'FontSize',fontsize, ...
              'HorizontalAlignment','right','VerticalAlignment','bottom');
  else
    i = find(mm == mm(end));
    ht = text(median(jd1(i)),yt,lmonth{mm(end)},'FontSize',fontsize, ...
              'HorizontalAlignment','center','VerticalAlignment','bottom');
  end
  hg = [hg;ht];
end


%-------------------------------------------------------------------------------

function hg = months(jd)
%MONTHS Label months.

global MTICK SEP MONTH

smonth = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct', ...
          'Nov','Dec'};
fontsize = get(gca,'FontSize');

[yy,mm,dd] = gregorian(jd);
yy = repmat(yy(1):yy(2),12,1);
mm = repmat([1:12]',1,size(yy,2));
yy = yy(:);
mm = mm(:);
xtick = julian(yy,mm,ones(size(mm)));
i = find(xtick >= jd(1) & xtick <= jd(2));
xtick = xtick(i);
yy = yy(i);
mm = mm(i);

ilabel = [1:4];
[dummy,i] = min(abs(diff(jd)/MONTH./ilabel - MTICK));
inc = ilabel(i);

for i = 1:length(xtick)
  if rem(mm(i)-1,inc) == 0
    xticklabel{i} = smonth{mm(i)};
  else
    xticklabel{i} = '';
  end   
end

set(gca,'XTick',xtick,'XTickLabel',xticklabel)

% label years
i = find(mm == 7 & xtick >= jd(1) & xtick <= jd(2));
xt = xtick(i);
yt = yloc(SEP,1);
hg = text(xt,yt*ones(size(xt)),int2str(yy(i)),'FontSize',fontsize, ...
          'HorizontalAlignment','center','VerticalAlignment','bottom');

% first year
[yy0,mm0,dd0] = gregorian(xtick(1));
if (mm0 > 7 & mm0 < 12-inc+1)
  if yy(end) - yy(1) > 1
    ht = text(jd(1),yt,int2str(yy(1)),'FontSize',fontsize, ...
	      'HorizontalAlignment','left','VerticalAlignment','bottom');
  else
    i = find(yy == yy(1));
    ht = text(median(xtick(i)),yt,int2str(yy(1)),'FontSize',fontsize, ...
	      'HorizontalAlignment','center','VerticalAlignment','bottom');
  end
  hg = [hg;ht];
end

% last year
[yy0,mm0,dd0] = gregorian(xtick(end));
if (mm0 < 7 & mm0 > inc)
  if yy(end) - yy(1) > 1
    ht = text(jd(2),yt,int2str(yy(end)),'FontSize',fontsize, ...
	      'HorizontalAlignment','right','VerticalAlignment','bottom');
  else
    i = find(yy == yy(end));
    ht = text(median(xtick(i)),yt,int2str(yy(end)),'FontSize',fontsize, ...
	      'HorizontalAlignment','center','VerticalAlignment','bottom');
  end
  hg = [hg;ht];
end

% bars
i = find(mm == 1 & xtick >= jd(1) & xtick <= jd(2));
xt = xtick(i);
yt = yloc(SEP,1);
l = length(xt);
if l > 0
  ht = text(xt,yt*ones(1,l),repmat('|',l,1),'FontSize',fontsize, ...
       'HorizontalAlignment','center','VerticalAlignment','bottom');
  hg = [hg;ht];
end


%-------------------------------------------------------------------------------

function years(jd)
%YEARS Label years.

global MTICK YEAR

ilabel = [1,2,5,10];
[dummy,i] = min(abs(diff(jd)/YEAR./ilabel - MTICK));
inc = ilabel(i);

[yy,mm,dd] = gregorian(jd);
yy = yy(1):yy(2);
xtick = julian(yy,ones(size(yy)),ones(size(yy)));

for i = 1:length(xtick)
  if rem(yy(i),inc) == 0
    xticklabel{i} = sprintf('%4d',yy(i));
  else
    xticklabel{i} = '';
  end   
end

set(gca,'XTick',xtick,'XTickLabel',xticklabel)


%-------------------------------------------------------------------------------

function decades(jd)
%DECADES Label decades.

global MTICK YEAR

ilabel = [5,10,20,50,100];
[dummy,i] = min(abs(diff(jd)/YEAR./ilabel - MTICK));
inc = ilabel(i);

[yy,mm,dd] = gregorian(jd);
yy = [round(yy(1)/inc) round(yy(2)/inc)]*inc;
yy = yy(1):inc:yy(end);
xtick = julian(yy,ones(size(yy)),ones(size(yy)));

for i = 1:length(xtick)
  xticklabel{i} = sprintf('%4d',yy(i));
end

set(gca,'XTick',xtick,'XTickLabel',xticklabel)


%-------------------------------------------------------------------------------

function y = yloc(dy,n)
%YLOC Determine y-location of labels.
%  Y = YLOC(DY) returns N y-locations of labels using the font height
%  separated by DY points.

n = n + 1;

% save current axes units
units = get(gca,'Units');

% y-position of x-axis in normalized units
set(gca,'Units','normalized')
ylim = get(gca,'YLim');

% position of axes and font size in points
set(gca,'Units','points')
pos = get(gca,'Position');
fontsize = get(gca,'FontSize');

% location of labels
ydir = get(gca,'ydir');
if strcmp(ydir,'reverse')
  y = ylim(2) + [1:n]*(fontsize + dy)*diff(ylim)/pos(4);
elseif strcmp(ydir,'normal')
  y = ylim(1) - [1:n]*(fontsize + dy)*diff(ylim)/pos(4);
end
% restore axes units
set(gca,'Units',units)

y = y(2:n);



