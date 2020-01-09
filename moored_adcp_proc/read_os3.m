function d = read_os3(file, varargin)
%function d = read_os(file, range, step, 'vel', ...)

% 00/01/31 EF
% 00/07/05 JH fixed ipp bug (here and in read_bb.m)
% 00/08/10 EF fixed JH's fix; added pg; fixed NaN output for missing data
% 00/09/20 EF Considerable additions and changes on Kaiyo, including
%             access to nav data type.
 
global id out yearbase

id_arg = 0;      % id is global so it only has to be called once.
                 % Call with 1 instead of zero to get NB if both NB
                 % and BB are being collected.
yearbase = 0;    % default: use first year in the file
ens_list = [];
step = 1;
start = -1;
 stop = 1000000000;
out.vel = 0; out.cor = 0; out.amp = 0; out.ts = 0; out.bt = 0; out.pg = 0;
out.nav = 0; out.ends = 0;

get_numeric = 0;
for iarg = 1:(nargin-1)
   arg = varargin{iarg};
   if ischar(arg)
      get_numeric = 0;
      switch lower(arg)
         case 'ends', out.ends = 1;
         case 'ens_list', get_numeric = 1;
         case 'yearbase', get_numeric = 2;
         case 'second_set', id_arg = 1;  % select NB when both are collected
         case 'vel', out.vel = 1;
         case 'cor', out.cor = 1;
         case 'amp', out.amp = 1;
         case 'pg', out.pg = 1;
         case 'ts', out.ts = 1;
         case 'bt', out.bt = 1;
         case 'nav', out.nav = 1; % not always present, so not in "all"
         case 'all', out.vel = 1; out.cor = 1; out.amp = 1;
                     out.ts = 1; out.bt = 1; out.pg = 1;
         otherwise
            error('invalid data type requested')
      end
   else   %numeric
      if get_numeric == 1,
         ens_list = arg;
         get_numeric = 0;
      elseif get_numeric == 2,
         yearbase = arg;
         get_numeric = 0;
      else
         nn = length(arg);
         if nn == 1,
            step = arg;
         elseif nn == 2,
            start = arg(1);
            stop = arg(2);
         else
            error('Invalid size of numeric argument')
         end
      end
   end
end

id = os_id(id_arg);

fid = fopen(file, 'r', 'ieee-le');
[config, buf] = read_buf(fid);  % First time: read config.
if ischar(buf)
   d.error = buf;
   return
end

fseek(fid, 0, 'eof');
len_file = ftell(fid);
nbytes = config.nbytes + 2;
n_profiles = floor(len_file / nbytes);
first_profile = buf.ens_num;
last_profile = first_profile + n_profiles - 1;


start = max([start, first_profile]);
stop = min([stop, last_profile]);

if out.ends & n_profiles > 4,
   ens_list = [first_profile, (first_profile + 1), (last_profile - 1), last_profile];
elseif isempty(ens_list)
   ens_list = start:step:stop;
end

NP = length(ens_list);  %%1 + floor((stop - start) / step);
NC = config.ncells;
NB = 4; % beams

if out.vel, d.vel = zeros(NC, NB, NP); end
if out.cor, d.cor = zeros(NC, NB, NP); end
if out.amp, d.amp = zeros(NC, NB, NP); end
if out.pg,  d.pg  = zeros(NC, NB, NP); end
if out.ts,
   d.heading = zeros(NP, 1);
   d.pitch =      d.heading;
   d.roll =       d.heading;
   d.soundspeed = d.heading;
   d.temperature =d.heading;
   d.pressure    =d.heading;
   d.juliandate  =d.heading;
   
end
if out.bt,
   d.bt.vel = zeros(NB, NP);
   d.bt.range = d.bt.vel;
   d.bt.cor =   d.bt.vel;
   d.bt.amp =   d.bt.vel;
   d.bt.rssi =  d.bt.vel;
end
if out.nav,
   d.nav.sec_pc_minus_utc = zeros(NP, 1);
   d.nav.txy1 = zeros(3, NP);
   d.nav.txy2 = d.nav.txy1;
end

% Always include dday and ens_num.
d.dday = zeros(NP, 1);
d.ens_num = zeros(NP, 1);
d.num_pings = zeros(NP, 1);
d.config = config;

d.depth = config.tr_depth + config.bin1distance + (0:(NC-1)) * config.cell;

ipp = 1;
for ip = ens_list
    
 %  try
   fseek(fid, nbytes * (ip - first_profile), 'bof');
   [config, buf] = read_buf(fid, config);
   if ischar(buf),
      d.error = buf;
      break
   end
   d.dday(ipp) = buf.dday;
   d.ens_num(ipp) = buf.ens_num;
   d.num_pings(ipp) = buf.num_pings;
   if out.ts,
      d.heading(ipp) = buf.heading;
      d.pitch(ipp)   = buf.pitch;
      d.roll(ipp)    = buf.roll;
      d.soundspeed(ipp)  = buf.soundspeed;
      d.temperature(ipp) = buf.temperature;
      d.pressure(ipp) = buf.pressure;
      d.juliandate(ipp) = buf.juliandate;
   end
   if out.vel, d.vel(:, :, ipp) = buf.vel; end
   if out.cor, d.cor(:, :, ipp) = buf.cor; end
   if out.amp, d.amp(:, :, ipp) = buf.amp; end
   if out.pg,  d.pg(:, :, ipp) = buf.pg; end
   if out.bt,
      d.bt.vel(:, ipp) = buf.bt.vel;
      d.bt.range(:, ipp) = buf.bt.range;
      d.bt.cor(:, ipp) = buf.bt.cor;
      d.bt.amp(:, ipp) = buf.bt.amp;
      d.bt.rssi(:, ipp) = buf.bt.rssi;
   end
   if out.nav,
      d.nav.sec_pc_minus_utc(ipp) = buf.nav.sec_pc_minus_utc;
      d.nav.txy1(:,ipp) = buf.nav.txy1;
      d.nav.txy2(:,ipp) = buf.nav.txy2;
   end
   ipp = ipp + 1;
   %end
end

if out.ts,
   d.heading = 0.01 * d.heading;
   d.pitch = 0.01 * to_signed(d.pitch);
   d.roll = 0.01 * to_signed(d.roll);
   d.temperature = 0.01 * to_signed(d.temperature);
   d.pressure = 0.001 * to_signed_long(d.pressure);
end

if out.vel, d.vel = 0.001 * to_signed(d.vel); end

if out.bt,
   d.bt.vel = 0.001 * to_signed(d.bt.vel);
   d.bt.range = 0.01 * d.bt.range;  % BAD is 0; leave it.
end


fclose(fid);
d.error = [];




%=========================================================================
function [config, buf] = read_buf(fid, config)

global id out yearbase

if nargin == 1, config.yearbase = yearbase; end % initializes config
start = ftell(fid);
[ss, count] = fread(fid, 2, 'uint16');
if count ~= 2,
   buf = sprintf('error 1, only %d bytes at start of buffer\n', count);
   return
end
if ss(1) ~= id.header,    % 7F7F
   buf = 'error 2, Header ID/Data Source not found in read_buf';
   return
end
config.nbytes = ss(2);  % Total bytes up to, not including, checksum.
fseek(fid, start, 'bof');
[rec, count] = fread(fid, config.nbytes, 'uint8');
if count ~= config.nbytes,
   buf = sprintf('error 3, only %d bytes read, %d expected\n', ...
                     count, config.nbytes);
   return
end
test_cs = mod(sum(rec), 65536);
cs = fread(fid, 1, 'uint16');
if test_cs ~= cs,
   buf = ['error 4, Checksum failure; CS is ', int2str(cs), ' versus calculated ', ...
            int2str(test_cs)];
   return
end

if nargin == 1,
   config.ndata = rec(6);
   ii = 7 + 2 * (0:(config.ndata-1));
   config.offsets = rec(ii) + 256 * rec(ii+1);
   config.ids = rec(config.offsets + 1) + rec(config.offsets + 2) * 256;
   %%% fixed leader:
   ofs = config.offsets(find(config.ids == id.flead));
   config.nbeams_used = rec(ofs + 9);
   config.nbeams_used = config.nbeams_used(1);
   config.ncells = rec(ofs + 10);
   config.ncells = config.ncells(1);
   config.cell = 0.01 * (rec(ofs + 13) + 256 * rec(ofs + 14));
   config.cell =config.cell(1);
   config.blank = 0.01 * (rec(ofs + 15) + 256 * rec(ofs + 16));
   config.blank=config.blank(1);
   config.mode = rec(ofs + 17);
   config.mode=config.mode(1);
   config.code_reps = rec(ofs + 19);
   config.code_reps=config.code_reps(1);
   config.EZ = rec(ofs + 31);
   config.EZ=config.EZ(1);
   config.bin1distance = 0.01 * (rec(ofs + 33) + 256 * rec(ofs + 34));
   config.bin1distance=config.bin1distance(1);
   config.pulse = 0.01 * (rec(ofs + 35) + 256 * rec(ofs + 36));
   config.pulse=config.pulse(1);
   config.lag = 0.01 * (rec(ofs + 41) + 256 * rec(ofs + 42));
   config.lag=config.lag(1);
   config.head_align = 0.01 * (rec(ofs + 27) + 256 * rec(ofs + 28));
   config.head_align=config.head_align(1);
   config.head_bias = 0.01 * (rec(ofs + 29) + 256 * rec(ofs + 30));
   config.head_bias=config.head_bias(1);
   lsb = rec(ofs + 5);
   msb = rec(ofs + 6);
   config.sysconfig = decode_sysconfig(lsb, msb);
end

%%% fixed leader, but possibly variable data:
ofs = config.offsets(find(config.ids == id.flead));
buf.num_pings = rec(ofs + 11) + 256 * rec(ofs + 12);
buf.num_pings=buf.num_pings(1);
%%% variable leader:
ofs = config.offsets(find(config.ids == id.vlead));
if nargin == 1,
   config.tr_depth = (rec(ofs + 17) + 256 * rec(ofs + 18))/10;
end
buf.ens_num = rec(ofs + 3) + 256 * rec(ofs + 4) + 65536 * rec(ofs + 12);

tt = rec(ofs + (5:11));
tt(1) = tt(1) + 1900;
if (tt(1) < 1980)  tt(1) = tt(1) + 100; end
 buf.juliandate = julian(tt(1),tt(2),tt(3),tt(4)+(tt(5)+tt(6)/100)/60);

if ~config.yearbase, config.yearbase = tt(1); end
try
buf.dday = julian(config.yearbase,tt(1),tt(2),tt(3),tt(4),tt(5),tt(6))+ tt(7) / 8640000;
catch
  buf.dday = NaN;
end

buf.ymdhmsh = tt;
if out.ts
   buf.heading = (rec(ofs + 19) + 256 * rec(ofs + 20));
   buf.pitch   = (rec(ofs + 21) + 256 * rec(ofs + 22));
   buf.roll    = (rec(ofs + 23) + 256 * rec(ofs + 24));
   buf.pressure = (rec(ofs + 49) + 256 * rec(ofs + 50) + 65536 * rec(ofs + 51) + 65536*256 * rec(ofs + 52));
   buf.temperature = rec(ofs + 27) + 256 * rec(ofs + 28);
   buf.soundspeed = rec(ofs + 15) + 256 * rec(ofs + 16);
%    buf.year = rec(ofs + 58)*100 + rec(ofs + 59);
%    buf.month = rec(ofs + 60);
%    buf.day = rec(ofs + 61);
%    buf.hour = rec(ofs+62)+(rec(ofs+63)/60 + rec(ofs+64)/3600 + rec(ofs+65)/360000);


end

% In the following blocks, the "min" should not be needed
% when finding the offsets; we may eliminate it.
% The empty check is needed, however.

%%% velocity:
if out.vel
   ofs = min(config.offsets(find(config.ids == id.vel)));
   ii = 2 + (1:2:(config.ncells*8));
   if isempty(ofs),
      a = NaN * ii;
   else
      a = rec(ofs + ii) + 256 * rec(ofs + 1 + ii);
   end
   buf.vel = reshape(a, 4, config.ncells).';
   % We are deferring the handling of negative numbers and NaNs.
end

%%% correlation:
if out.cor
   ofs = min(config.offsets(find(config.ids == id.cor)));
   ii = 2 + (1:(config.ncells*4));
   if isempty(ofs)
      buf.cor = reshape(NaN * ii, 4, config.ncells).';
   else
      buf.cor = reshape(rec(ofs + ii), 4, config.ncells).';
   end
end

%%% amp:
if out.amp
   ofs = min(config.offsets(find(config.ids == id.amp)));
   ii = 2 + (1:(config.ncells*4));
   if isempty(ofs)
      buf.amp = reshape(NaN * ii, 4, config.ncells).';
   else
      buf.amp = reshape(rec(ofs + ii), 4, config.ncells).';
   end
end

%%% pg:
if out.pg
   ofs = min(config.offsets(find(config.ids == id.pg)));
   ii = 2 + (1:(config.ncells*4));
   if isempty(ofs)
      buf.pg = reshape(NaN * ii, 4, config.ncells).';
   else
      buf.pg = reshape(rec(ofs + ii), 4, config.ncells).';
   end
end

%% BT:
if out.bt
   ofs = min(config.offsets(find(config.ids == id.bt)));
   if isempty(ofs)

      config.bt.BP = NaN;
      config.bt.BC = NaN;
      config.bt.BA = NaN;
      config.bt.BX = NaN;

      ii = (0:2:6);
      buf.bt.range = NaN*ii;
      buf.bt.vel =   NaN*ii;
      buf.bt.cor =   NaN*ii;
      buf.bt.amp =   NaN*ii;
      buf.bt.rssi =  NaN*ii;

   else
      config.bt.BP = rec(ofs + 3) + 256 * rec(ofs + 4);
      config.bt.BC = rec(ofs + 7);
      config.bt.BA = rec(ofs + 8);
      config.bt.BX = rec(ofs + 71) + 256 * rec(ofs + 72);

      ii = (0:2:6); iii = 0:3;
      buf.bt.range = rec(ofs + ii + 17) + 256 * rec(ofs + ii + 18) + ...
                     65536 * rec(ofs + iii + 78);

      buf.bt.vel = rec(ofs + ii + 25) + 256 * rec(ofs + ii + 26);
      buf.bt.cor = rec(ofs + iii + 33);
      buf.bt.amp = rec(ofs + iii + 37);
      buf.bt.rssi = rec(ofs + iii + 73);
   end
end

%% nav:
if out.nav
   ofs = min(config.offsets(find(config.ids == id.nav)));
   if isempty(ofs)
      buf.nav.sec_pc_minus_utc = NaN;
      buf.nav.txy1 = NaN * ones(3,1);
      buf.nav.txy2 = buf.nav.txy1;
   else
      buf.nav.sec_pc_minus_utc = 1e-3 * buf_long(rec(ofs + (11:14)));


      day_base = to_day(config.yearbase, ...
                        rec(ofs + 5) + 256 * rec(ofs + 6), ...
                        rec(ofs + 4), rec(ofs + 3));
      t1d = day_base + (1e-4 / 86400) * buf_long(rec(ofs + (7:10)));
      t2d = day_base + (1e-4 / 86400) * buf_long(rec(ofs + (23:26)));

      buf.nav.txy1(1) = t1d + round(buf.dday - t1d);
      buf.nav.txy2(1) = t2d + round(buf.dday - t2d);
      buf.nav.txy1(2:3) = decode_fix(rec(ofs + (15:22)));
      buf.nav.txy2(2:3) = decode_fix(rec(ofs + (27:34)));
   end
end

%--------------------------------------------------------
function c = decode_sysconfig(l, m)

freqs = [75, 150, 300, 600, 1200, 2400, 38];
freqbits = bitand(l, 7);
c.frequency = freqs(1 + freqbits);
c.up =  bitand(bitshift(l, -7), 1);
c.convex = bitand(bitshift(l, -3), 1);

angles = [15, 20, 30, NaN];
anglebits = bitand(m, 3);
c.angle = angles(1 + anglebits);



%---------------------------------------------
function s = to_signed(u)
s = u;
ineg = u > 32767;
s(ineg) = s(ineg) - 65536;
ibad = s == -32768;
s(ibad) = NaN;

%---------------------------------------------
function s = to_signed_long(u)
s = u;
ineg = u > 2147483647 ;
s(ineg) = s(ineg) - 4294967296;
ibad = s == -2147483648;
s(ibad) = NaN;

%---------------------------------------------
function s = buf_long(a)
s = a(1) + 256 * a(2) + 65536 * a(3) + 16777216 * a(4);
s = to_signed_long(s);

%---------------------------------------------
function xy = decode_fix(r)
xy = [buf_long(r(5:8)); buf_long(r(1:4))] * 180 / 2147483648;



%----------------------------------------------
function id = os_id(incr)
id.header = hex2dec('7f7f');
id.flead = 0 + incr;
id.vlead = hex2dec('0080') + incr;
id.vel = hex2dec('0100') + incr;
id.cor = hex2dec('0200') + incr;
id.amp = hex2dec('0300') + incr;
id.pg = hex2dec('0400') + incr;
id.status = hex2dec('0500') + incr; % Useless; ignore it.
id.bt = hex2dec('0600') + incr;
% Added in STA and LTA files:
id.nav = hex2dec('2000');           % No separate nav type.
id.fix_att = hex2dec('3000');       % Or fixed attitude.
id.var_att1 = hex2dec('30e0');
id.var_att2 = hex2dec('30e8');
id.var_att = id.var_att1; % default: NP only
                  % other possibilities to be added.














