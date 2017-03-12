function calmat = F_calmat(filename,config)
% This matlab function reads interferogram from OPUS file (filename), does
% DC correction, filters out bad ifg, and performs fft to derive spectrum
% using different apodization functions specified by config.apokind:
%          1: boxcar
%          2: triag
%          3: Happ-Genzel
%          4/5: Blackmann-Harris 3-term/4-term
%          6/7/8: Norton-Beer weak/medium/strong
% This function is inspired by Calpy Python package by Matthaus Kiel
% The matlab OPUS reading routine is kindly provided by Jacob Filik
% Written by Kang Sun on 2017/01/04

if ~isfield(config,'apokind')
    apokind = 7;
else
    apokind = config.apokind;
end

if ~isfield(config,'localt')
    localt = 0;
else
    localt = config.localt;
end

if ~isfield(config,'loclat')
    loclat = -999;loclon = -999;
else
    loclat = config.loclat;loclon = config.loclon;
end

if ~isfield(config,'locname')
    locname = 'no idea!!';
else
    locname = config.locname;
end

if ~isfield(config,'timecorr')
    timecorr = 0;
else
    timecorr = config.timecorr;
end

if ~isfield(config,'timezone')
    timezone = -5;
else
    timezone = config.timezone;
end

if ~isfield(config,'fieldstopdiam')
    fieldstopdiam = 0.6;
else
    fieldstopdiam = config.fieldstopdiam;
end

if ~isfield(config,'focallength')
    focallength = 127;
else
    focallength = config.focallength;
end

if ~isfield(config,'lowerlim')
    lowerlim = 0;upperlim = inf;
else
    lowerlim = config.lowerlim;upperlim = config.upperlim;
end

phasres = 0;   % choose the smallest possible phase resolution
ssords = 1;    % double-sided ifg
timingtype = 0;% if 0, timing at the start of the scan; if 1, at the end
calmat = [];   % output of this function

if apokind < 0
    try
        [spec,x,paramas] = ImportOpus(filename,'SampleSpectrum');
        warning('off','MATLAB:pfileOlderThanMfile');
        fclose('all');
    catch
        disp('Error: something went wrong when reading spectrum!):')
        disp(['Skipping file: ',filename])
        return
    end
    calmat.spec = spec(x >= lowerlim & x <= upperlim);
    calmat.x = x(x >= lowerlim & x <= upperlim);
else
    [ifg,~,paramas] = ImportOpus(filename,'SampleInterferogram');
    warning('off','MATLAB:pfileOlderThanMfile');
    fclose('all');
    
    if rem(paramas.Instrument.GFW+paramas.Instrument.GBW,2) ~= 0
        disp('Error: Number of good forward and backward scans is odd!):')
        disp(['Skipping interferogram: ',filename])
        return
    end
    
    nifg = length(ifg);
    ifg(nifg/2+1:end) = flipud(ifg(nifg/2+1:end));
    
    if sum(ifg > 0) > 0 && sum(ifg < 0) > 0
        disp('Values change sign. Really DC-Interferogram?')
        disp(['Skipping interferogram: ',filename])
        return
    end
    if sum(ifg == 0) > 0
        disp('Zero-entries in IFG. Really DC-Interferogram?')
        disp(['Skipping interferogram: ',filename])
        return
    end
    if rem(nifg,2) ~= 0
        disp('Error: number of entries in IFG is odd (damaged IFG?!):')
        disp(['Skipping interferogram: ',filename])
        return
    end
    disp('-> plausibility check passed!')
    
    % DC correction:
    % find the smoothed DC baseline
    ifgLeft = double(ifg(1:nifg/2));
    ifgRight = double(ifg(nifg/2+1:end));
    ifgLeft_smooth = F_ifg_smooth(ifgLeft);
    ifgRight_smooth = F_ifg_smooth(ifgRight);
    
    ifg_smooth = [ifgLeft_smooth(:);ifgRight_smooth(:)];
    calmat.ifg_smooth = ifg_smooth;
    
    % find ZPD positions and other stuff
    [maxAbsLeft, posMaxAbsLeft] = max(abs(ifgLeft));
    [minAbsLeft, posMinAbsLeft] = min(abs(ifgLeft));
    
    [maxAbsRight, posMaxAbsRight] = max(abs(ifgRight));
    [minAbsRight, posMinAbsRight] = min(abs(ifgRight));
    
    [maxAbsAroundLeftZPD, posMaxAbsAroundLeftZPD]...
        = max(abs(ifgLeft(posMaxAbsLeft-20:posMaxAbsLeft + 20)));
    [minAbsAroundLeftZPD, posMinAbsAroundLeftZPD]...
        = min(abs(ifgLeft(posMaxAbsLeft-20:posMaxAbsLeft + 20)));
    
    [maxAroundLeftZPD, posMaxAroundLeftZPD]...
        = max(ifgLeft(posMaxAbsLeft-20:posMaxAbsLeft + 20));
    [minAroundLeftZPD, posMinAroundLeftZPD]...
        = min(ifgLeft(posMaxAbsLeft-20:posMaxAbsLeft + 20));
    
    [maxAbsAroundRightZPD, posMaxAbsAroundRightZPD]...
        = max(abs(ifgRight(posMaxAbsRight-20:posMaxAbsRight + 20)));
    [minAbsAroundRightZPD, posMinAbsAroundRightZPD]...
        = min(abs(ifgRight(posMaxAbsRight-20:posMaxAbsRight + 20)));
    
    [maxAroundRightZPD, posMaxAroundRightZPD]...
        = max(ifgRight(posMaxAbsRight-20:posMaxAbsRight + 20));
    [minAroundRightZPD, posMinAroundRightZPD]...
        = min(ifgRight(posMaxAbsRight-20:posMaxAbsRight + 20));
    
    % calculation of mean values around ZPDs
    pointsArroundZPD = 100;
    ifgTemp = ifg_smooth;
    
    meanWhole = mean(ifgTemp);
    maxWhole = max(ifgTemp);
    minWhole = min(ifgTemp);
    percentDeviation = ((maxWhole - minWhole) / meanWhole) * 100.0;
    ifgTempLeft = ifgTemp((posMaxAbsLeft - pointsArroundZPD):...
        (posMaxAbsLeft + pointsArroundZPD));
    meanLeftZPD = mean(ifgTempLeft);
    ifgTempRight = ifgTemp((posMaxAbsRight - pointsArroundZPD):...
        (posMaxAbsRight + pointsArroundZPD));
    meanRightZPD = mean(ifgTempRight);
    calmat.ZPDRight = ifgRight((posMaxAbsRight - pointsArroundZPD):...
        (posMaxAbsRight + pointsArroundZPD));
    calmat.ZPDLeft = ifgLeft((posMaxAbsLeft - pointsArroundZPD):...
        (posMaxAbsLeft + pointsArroundZPD));
    calmat.percentDeviation = percentDeviation;
    
    % calculate ifg quality
    meanZPD = 0.5 * (meanLeftZPD + meanRightZPD);
    sizeIFGPeakLeft = maxAbsLeft - minAbsAroundLeftZPD;
    sizeIFGPeakRight = maxAbsRight - minAbsAroundRightZPD;
    s = sum(abs(log(abs(ifgTemp / meanZPD)))/length(ifgTemp));
    ifgQuality = s / (sizeIFGPeakLeft * sizeIFGPeakRight);
    if (abs(posMaxAbsRight - posMaxAbsLeft) > 20)
        % Max-Position is not symmetric -> bad IFG
        ifgQuality = - 2;
    end
    disp(['-> calculation of IFG-quality parameter: ',num2str(ifgQuality)])
    calmat.ifgQuality = ifgQuality;
    
    % perform DC correction
    ifgC = (ifg * meanZPD) ./ ifg_smooth - meanZPD;
    calmat.ifgC = ifgC;
    
    % determine position of centerburst
    nifgscan = nifg / 2;
    [~,nburstfwd] = max(abs(ifgC(1:nifgscan)));
    [~,nburstbwd] = max(abs(ifgC(nifgscan+1:end)));
    % matlab is 1 based, minus 1 for number of channels before burst
    nburstfwd = nburstfwd-1;
    nburstbwd = nburstbwd-1;
    % disp([ '   centerburst location (fwd): ', num2str(nburstfwd)]);
    % disp([ '   centerburst location (bwd):', num2str(nburstbwd)]);
    
    nradiusfwd = min(nburstfwd, nifgscan - nburstfwd - 1);
    nradiusbwd = min(nburstbwd, nifgscan - nburstbwd - 1);
    nradius = min(nradiusfwd, nradiusbwd);
    
    % disp([ '   nifgscan:', num2str(nifgscan)])
    % disp([ '   nburstfwd:', num2str(nburstfwd)])
    % disp([ '   nburstbwd:', num2str(nburstbwd)])
    nss = min(nifgscan - nburstfwd - 1,nifgscan - nburstbwd - 1);
    calmat.nburstfwd = nburstfwd;
    calmat.nburstbwd = nburstbwd;
    calmat.nradius = nradius;
    calmat.nss = nss;
    calmat.nifgscan = nifgscan;
    
    % preparing FFT
    nuesampling = paramas.Instrument.HFL;
    if paramas.Instrument.LFL > 0
        warning(['Warning!!! Low frequency limit > 0 but 0 is expected -->',...
            'wrong results!']);
    end
    
    % starting FWD FFT
    [specrepcfwd, phasrefwd, phasimfwd,nmax] =...
        F_fftmain(ifgC(1:length(ifgC)/2), nuesampling, nburstfwd,...
        apokind, ssords, nss, nradius, phasres);
    
    % starting BWD FFT
    [specrepcbwd, phasrebwd, phasimbwd,nmax] =...
        F_fftmain(ifgC(length(ifgC)/2+1:length(ifgC)), nuesampling,...
        nburstfwd, apokind, ssords, nss, nradius, phasres);
    
    spec = 0.5 * specrepcfwd + 0.5 * specrepcbwd;
    x = linspace(0,nuesampling,nmax/2+1);
    x(1) = 10.0e-10; % to avoid dividing through zero
    x = x(1:nmax/2 + 1);
    x = x(:);
    calmat.nuesampling = nuesampling;
    calmat.nmax = nmax;
    calmat.x = x(x >= lowerlim & x <= upperlim);
    calmat.spec = spec(x >= lowerlim & x <= upperlim);
end
% read time, calculate solar zenith angle and air mass
DAT = sscanf(paramas.SampleDataInterferogram.DAT(1:end-2),'%f/%f/%f');
TIM = sscanf(paramas.SampleDataInterferogram.TIM(1:12),'%f:%f:%f');
DAT = DAT(:)';TIM = TIM(:)';
tutc = datenum([DAT(end:-1:1),TIM])+timecorr/24;
tloc = tutc+timezone/24;
halfduration = paramas.Instrument.DUR/2/86400;
if timingtype == 0
    tutc = tutc+halfduration;
    tloc = tloc+halfduration;
elseif timingtype == 1
    tutc = tutc-halfduration;
    tloc = tloc+halfduration;
end
calmat.tutc = tutc;
calmat.tloc = tloc;
if ~isfield(config,'refraction')
    disp('Tell me you want refraction or not!')
    return % done
end
if ~isfield(config,'pathpt')
    config.refraction = 0;
    config.pathpt = [];
end
dtv = datevec(tutc); % stands for "date time vector"
pos = [loclon,loclat,localt];

[error, elev_unref, elev_ref, azim, am] = ...
    F_sun_position(dtv, pos, config.refraction,config.pathpt, 'SUN');
if error
    disp('An error occurred while calculating solar position.')
    disp('Check sun_pos2_log.txt and system compatibility')
else
    calmat.elev_unref = elev_unref;
    calmat.elev_ref = elev_ref;
    calmat.azim = azim;
    calmat.am = am;
end

function [error, elev_unref, elev_ref, azim, am] = ...
    F_sun_position(dtv,pos,refraction,pathpt,sunormoon)
fid = fopen('sun_pos_2.inp', 'w');
fprintf(fid,'%s',sprintf(['Input-File for Sun_Pos_2.exe\n', '\n', ...
    'Note: The Earth-Orientation parameter file can be downloaded here:',...
    '\n', ...
    'http://hpiers.obspm.fr/iers/eop/eopc04/eopc04.62-now', '\n',...
    '\n', ...
    'Needed input:', '\n', ...
    '- Year (e.g. 2010) (if only 2 digits, then <80 is interpreted as 20xx, >=80 as 19xx)', '\n',...
    '- Month (1-12)', '\n', ...
    '- Day', '\n', ...
    '- Hour (0-23)[Integer] (can also be from -1 to 47)', '\n', ...
    '- Minute (0-59)[Integer]', '\n', ...
    '- Second (0.0-60.0)[Float]', '\n', ...
    'Observer Position:', '\n', ...
    'Option 1:', '\n', ...
    '- latitude (degree, North is positive, South is negative, [Karlsruhe: 49.1])', '\n', ...
    '- longitude (degree, East is positive, West is negative, [Karlsruhe: 8.44])', '\n', ...
    '- altitude [km]', '\n', ...
    'Option 2:', '\n', ...
    '- name of station as given in station_list.inp. (E.g. Karlsruhe)', '\n\n', ...
    '$1', '\n', ...
    num2str(dtv(1),'%d'), '\n',...
    num2str(dtv(2),'%d'), '\n',...
    num2str(dtv(3),'%d'), '\n',...
    num2str(dtv(4),'%f'), '\n',...
    num2str(dtv(5),'%f'), '\n',...
    num2str(dtv(6),'%f'), '\n',...
    num2str(pos(2),'%f'), '\n',...
    num2str(pos(1),'%f'), '\n',...
    num2str(pos(3),'%f'), '\n\n\n\n', ...
    'Optional input:\n', ...
    '- Include refraction: (0:false, 1:true)\n',...
    '- path to "pt.prf"-file  Use . if it is in same folder than sun_pos_2.exe)\n',...
    '- SUN or MOON\n\n', ...
    '$2\n', num2str(refraction,'%d'), '\n', pathpt, '\n', sunormoon]));
fclose(fid);
error = system('sun_pos_2.exe > sun_pos2_log.txt 2>&1');
if error ~= 0
    elev_unref = 0;
    elev_ref = 0;
    azim = 0;
    am = 0;
    return
else
    fid = fopen('sun_pos_2.out','r');
    error = str2double(fgets(fid));
    elev_unref = str2double(fgets(fid));
    elev_ref = str2double(fgets(fid));
    azim = str2double(fgets(fid));
    fgets(fid);fgets(fid);fgets(fid);
    am = str2double(fgets(fid));
    fclose(fid);
end

function [specrepc, phasre, phasim,nmax] = ...
    F_fftmain(ifgacfwd,nuesampling,nburstfwd,apokind,ssords,nss,nradius,phasres)
% main function to prepare and call the fft function.
% ifgacfwd = ifgC(1:nifgscan);
if phasres > 0
    nradiusphas = int(2 * nuesampling / phasres);
    if nradiusphas > nradius
        nradiusphas = nradius;
        warning( 'Warning: Not enough pts for requested phase res!')
    end
else
    nradiusphas = nradius;
end

nmax = 2^nextpow2(nss+nradius);
%     while nmax < nss + nradius
%         nmax = nmax * 2;
%     end
if ssords == 0
    % single-sided
    %     disp([ 'OPDmax single-sided IFG:', num2str(nss / (2 * nuesampling))])
    if nradius > nss - 1
        nradius = nss - 1;
    end
    %     disp([ 'OPD available for phase calc:', num2str(nradius / (2 * nuesampling))])
    %     disp(['OPD requested for phase calc:', num2str(nradiusphas / (2 * nuesampling))])
    
    nmax = nmax * 2; % Apply additional zerofill of factor 2 for SS-IFG!!
    %     disp( ['FFT size: ', num2str(nmax)])
    deltanue = nuesampling / single(nmax / 2);
elseif ssords == 1
    % double-sided
    %     disp([ '   OPDmax: ',num2str((nradius / (2 * nuesampling)))])
    %     disp([ '   OPD requested for phase calc: ',...
    %         num2str((nradiusphas / (2 * nuesampling)))])
    %     disp( ['   FFT size: ', num2str(nmax)])
    deltanue = nuesampling / single(nmax / 2);
end
% disp('   initializing FFT arrays ...')
ifgre = zeros(nmax,1,'single');
%     ifgim = zeros(nmax,1,'single');
ifgphasre = zeros(nmax,1,'single');
%     ifgphasim = zeros(nmax,1,'single');

if ssords == 0
    %     disp( '   before FFT, preparing SS-IFG')
    xx = nburstfwd+1:nburstfwd + nss + 1;
    xx = xx(:);
    ifgre(1:length(xx)) = ifgacfwd(xx).* ...
        F_apoifg(apokind,nburstfwd,nradius,xx).*...
        F_aporamp(nburstfwd,nradiusphas,xx);
    ifgphasre(1:length(xx)) = ifgacfwd(xx).*...
        F_apophas(nburstfwd,nradiusphas,xx);
    
    xx = nburstfwd:-1:nburstfwd-nradiusphas+1;
    xx = xx(:);
    ifgre(end:-1:end-length(xx)+1) = ifgacfwd(xx).*...
        F_apoifg(apokind,nburstfwd,nradius,xx).*...
        F_aporamp(nburstfwd,nradiusphas,xx);
    ifgphasre(end:-1:end-length(xx)+1) = ifgacfwd(xx).*...
        F_apophas(nburstfwd,nradiusphas,xx);
    
elseif ssords == 1
    %     disp( '   before FFT, preparing DS-IFG')
    xx = nburstfwd+1:nburstfwd + nradius+1;
    xx = xx(:);
    ifgre(1:length(xx)) = 0.5*ifgacfwd(xx) .* ...
        F_apoifg(apokind,nburstfwd,nradius,xx);
    ifgphasre(1:length(xx)) = ifgacfwd(xx).*...
        F_apophas(nburstfwd,nradiusphas,xx);
    
    xx = nburstfwd:-1:nburstfwd-nradius+1;
    xx = xx(:);
    ifgre(end:-1:end-length(xx)+1) = 0.5*ifgacfwd(xx).*...
        F_apoifg(apokind,nburstfwd,nradius,xx);
    ifgphasre(end:-1:end-length(xx)+1) = ifgacfwd(xx).*...
        F_apophas(nburstfwd,nradiusphas,xx);
end
Y = fft(ifgphasre, nmax);
phasre = real(Y);
phasim = imag(Y);

Y = fft(ifgre, nmax);
specre = real(Y);
specim = imag(Y);

%     corrupt_array = 1.0./sqrt(phasre .* phasre + phasim .* phasim);
norm = 1.0 ./ sqrt(phasre .* phasre + phasim .* phasim);
phasre = phasre .* norm;
phasim = phasim .* norm;
%     specrepc  = zeros(nmax,1,'single');
specrepc = specre .* phasre + specim .* phasim;
specrepc = specrepc(1:nmax/2 + 1);
phasre = phasre(1:nmax/2 + 1);

function ifg_fwd_smooth = F_ifg_smooth(ifg_fwd)
% a neater way to calculate smooth ifg DC signal. avSize and smooth
% window can be further tuned.
% Written by Kang Sun on 2017/1/1

xx = 1:length(ifg_fwd);
avSize = 1000;
Edges = linspace(1,length(ifg_fwd),round(length(ifg_fwd)/avSize));
binmid = Edges(1:end-1)+(Edges(2)-Edges(1))/2;
[~,ind] = histc(xx,Edges);
yy1 = accumarray(ind(:),ifg_fwd,[],@mean);
yy1 = yy1(1:end-1);
binmid = [1;binmid(:);length(ifg_fwd)];
yy1 = [mean(ifg_fwd(1:200));yy1(:);mean(ifg_fwd(end-200:end))];
if exist('smooth','file')
    yy1_s = smooth(yy1,0.05);
else
    yy1_s = masmooth(1:length(yy1),yy1,round(length(yy1)*0.05),'mean',5);
end
ifg_fwd_smooth = spline(binmid,yy1_s,xx);

function apoifg = F_apoifg(apokind,nburst,nopdmax,xx)
xwert = abs(xx - nburst) / single(nopdmax);

if apokind == 1
    apoifg = ones(size(xx));
elseif apokind == 2
    apoifg = 1.0 - xwert;
elseif apokind == 3
    apoifg = 0.53856 + 0.46144 * cos(pi * xwert);
elseif apokind == 4
    apoifg = 0.42323 + 0.49755 * cos(pi * xwert) + 0.07922 * cos(2.0 * pi * xwert);
elseif apokind == 5
    apoifg =  0.35875 + 0.48829 * cos(pi * xwert) + 0.14128 * cos(2.0 * pi * xwert) + 0.01168 * cos(3.0 * pi * xwert);
elseif apokind == 6
    term = 1.0 - xwert .* xwert;
    apoifg =  0.384093 - 0.087577 * term + 0.703484 * term .* term;
elseif apokind == 7
    term = 1.0 - xwert .* xwert;
    apoifg = 0.152442 - 0.136176 * term + 0.983734 * term .* term;
elseif apokind == 8
    term = (1.0 - xwert .* xwert)  .* (1.0 - xwert .* xwert);
    apoifg =  0.045335 + 0.554883 * term + 0.399782 * term .* term;
end
apoifg(xx-nburst > nopdmax) = 0;

function apophas = F_apophas(nburst, ncutoff, xx)
apophas = zeros(size(xx));
temp = 0.5 + 0.5 * cos(pi * single(abs(xx - nburst)) / single(ncutoff));
apophas(abs(xx-nburst) < ncutoff) = temp(abs(xx-nburst) < ncutoff);

function wert = F_aporamp(nburst, nradius, xx)
wert = 0.5 + 0.5 * single(xx - nburst) / single(nradius);
wert(wert > 1) = 1;
wert(wert < 0) = 0;

