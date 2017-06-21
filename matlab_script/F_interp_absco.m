function Vout2 = F_interp_absco(inp)
% Function to interpolate absco table with dimension P, T, Broadener,
% Wavenumber into query points Pq, Tq, (Broadener not implemented), and Wq
% Previous failure attempts can be found in test_interpolation.m

% The absco table has to be on semi-aligned temperature grid, e.g.,
% [230,240,250,260,270,280,290,300;...
%      240,250,260,270,280,290,300,310]

% Written by Kang Sun on 2017/04/29

tempgrid = inp.tempgrid;
presgrid = inp.presgrid;
wavegrid = inp.wavegrid;
if isfield(inp,'bgrid')
    bgrid = inp.bgrid;
else
    bgrid = 0;
end

V = inp.absco;

% convert presgrid from Pa or atm to hPa if necessary
if max(presgrid) > 7.7e4 % I guess a good pressure profile extends below 770 hPa
    presgrid = presgrid/100;
end
if presgrid < 2
    presgrid = presgrid*1.01325e3;
end

Tq = inp.Tq;% query temperature in K
Pq = inp.Pq;% query pressure in hPa

if Pq <= min(presgrid) && Pq >= max(presgrid)
    error('Pressure out of range!')
end

if isfield(inp,'Wq')
    Wq = inp.Wq;
elseif isfield(inp,'wave_downsample_step')
    Wq = wavegrid(1:inp.wave_downsample_step:end);
else
    Wq = wavegrid;
end
Wq = Wq(:);

if isfield(inp,'Bq')
    Bq = inp.Bq;
else
    Bq = 0;
end

if isfield(inp,'n_interp_level')
    n_interp_level = inp.n_interp_level;
else
    n_interp_level = 3;
end

if isfield(inp,'T_ext')
    T_ext = inp.T_ext;
else
    T_ext = 20;
end
% trim temperature outliers
[~, I] = sort(abs(presgrid-Pq));
int_P = sort(I(1:n_interp_level));

% deal with the bug that pressures all on one side of Pq
if all(presgrid(int_P)-Pq > 0)
    I(n_interp_level) = I(1)-1;
    int_P = sort(I(1:n_interp_level));
end
if all(presgrid(int_P)-Pq < 0)
    I(n_interp_level) = I(1)+1;
    int_P = sort(I(1:n_interp_level));
end

if Tq <= max(tempgrid(I(1:n_interp_level),1))
    Tq0 = Tq;
    Tq = max(tempgrid(I(1:n_interp_level),1))+1;
    warning(['Temperature too low at ',num2str(Pq),' hPa! Changed from ',num2str(Tq0),' to ',num2str(Tq)])
end

if  Tq >= min(tempgrid(I(1:n_interp_level),end))
    Tq0 = Tq;
    Tq = min(tempgrid(I(1:n_interp_level),end))-1;
    warning(['Temperature too high at ',num2str(Pq),' hPa! Changed from ',num2str(Tq0),' to ',num2str(Tq)])
end

% make the array size more tangible
presgrid = presgrid(int_P);
presgrid = presgrid(:);
tempgrid = tempgrid(int_P,:);
V = V(:,:,:,int_P);
V = permute(V,[4 3 2 1]);
% let's make the dimension order clear: Pressure, temperature, water
% vapor, wavelength

T_ext_left = min(Tq-min(tempgrid,[],2));
T_ext_right = min(max(tempgrid,[],2)-Tq);

if T_ext_left <= 0 || T_ext_right <= 0
    error('Temperature out of range because pressure ext is too large!')
end

int_temp = false(size(tempgrid));
for ip = 1:length(presgrid)
    if T_ext_left < T_ext % left corner case
        int_temp(ip,:) = tempgrid(ip,:) >= Tq-T_ext_left & tempgrid(ip,:) <= Tq+T_ext;
    elseif T_ext_right < T_ext % right corner case
        int_temp(ip,:) = tempgrid(ip,:) >= Tq-T_ext & tempgrid(ip,:) <= Tq+T_ext_right;
    else
        int_temp(ip,:) = tempgrid(ip,:) >= Tq-T_ext & tempgrid(ip,:) <= Tq+T_ext;
    end
end

ntemp = min(sum(int_temp,2));
if ntemp ~= max(sum(int_temp,2));
    warning('Temperature interval should be equal for each level. You have a problem!')
end

tempgrid_trim = nan(length(presgrid),ntemp,'single');
V_trim = nan(size(V,1),ntemp,size(V,3),size(V,4),'single');

for ip = 1:length(presgrid)
    tempgrid_trim(ip,:) = tempgrid(ip,int_temp(ip,:));
    V_trim(ip,:,:,:) = V(ip,int_temp(ip,:),:,:);
end

% regular grid method, interpn
Vout2 = nan(length(bgrid),length(Wq),'single');
for iwv = 1:length(bgrid)
    Vtemp = squeeze(V_trim(:,:,iwv,:));
    Pg = repmat(presgrid,[1,ntemp,length(wavegrid)]);
    Tg = repmat(tempgrid_trim,[1 1 length(wavegrid)]);
    Wg = repmat(wavegrid,[1 length(presgrid) ntemp]);
    Wg = permute(Wg,[2 3 1]);
    Vout2(iwv,:) = interpn(single(Pg),single(Tg),single(Wg),single(Vtemp),...
        single(ones(size(Wq))*Pq),single(ones(size(Wq))*Tq),single(Wq));
end