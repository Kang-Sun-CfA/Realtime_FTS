function output = F_calmat_ils(input)
if isfield(input,'filename')
    filename = input.filename;
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
%     disp('-> plausibility check passed!')
    
    % DC correction:
    % find the smoothed DC baseline
    ifgLeft = double(ifg(1:nifg/2));
    ifgRight = double(ifg(nifg/2+1:end));
    ifgLeft_smooth = F_ifg_smooth(ifgLeft);
    ifgRight_smooth = F_ifg_smooth(ifgRight);
    
    ifg_smooth = [ifgLeft_smooth(:);ifgRight_smooth(:)];
    output.ifg_smooth = ifg_smooth;
    
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
    output.ZPDRight = ifgRight((posMaxAbsRight - pointsArroundZPD):...
        (posMaxAbsRight + pointsArroundZPD));
    output.ZPDLeft = ifgLeft((posMaxAbsLeft - pointsArroundZPD):...
        (posMaxAbsLeft + pointsArroundZPD));
    output.percentDeviation = percentDeviation;
    
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
%     disp(['-> calculation of IFG-quality parameter: ',num2str(ifgQuality)])
    output.ifgQuality = ifgQuality;
    
    % perform DC correction
    ifgC = (ifg * meanZPD) ./ ifg_smooth - meanZPD;
    
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
    
    % preparing FFT
    nuesampling = paramas.Instrument.HFL;
    if paramas.Instrument.LFL > 0
        warning(['Warning!!! Low frequency limit > 0 but 0 is expected -->',...
            'wrong results!']);
    end
    nuesampling = double(nuesampling);
    nradius = double(nradius);
else
    nuesampling = double(input.nuesampling);
    nradius = double(input.nradius);
end
    tspace = (-nradius:nradius)/2/nuesampling;
     sim_wn = input.sim_wn;
    phase_error = input.phase_error;
    zerofillingfactor = input.zerofillingfactor;
    ispace = cos(2*pi*sim_wn*(tspace)-phase_error);
    ispace = ispace(:);
%     ispace(tspace < 0) = flipud(ispace(tspace > 0));

    apokind = input.apokind;
%     plot(tspace,ispace)
%     xlim([-5/nuesampling,5/nuesampling])
    [specils, phasrefwd, phasimfwd,nmax,specre,specim] =...
        F_fftmain(ispace, nuesampling, nradius,...
        apokind, 1, nradius, nradius, 0,zerofillingfactor);
    x = linspace(0,nuesampling,nmax/2+1);
    x(1) = 10.0e-10; % to avoid dividing through zero
    x = x(1:nmax/2 + 1);
    x = x(:);
    output.ilsy = specils;
    output.ilsx = x;
    output.nmax = nmax;
    output.specre = specre;
    output.specim = specim;
function [specrepc, phasre, phasim,nmax,specre,specim] = ...
    F_fftmain(ifgacfwd,nuesampling,nburstfwd,apokind,ssords,nss,nradius,phasres,zerofillingfactor)
% main function to prepare and call the fft function.
% ifgacfwd = ifgC(1:nifgscan);
if ~exist('zerofillingfactor','var')
    zerofillingfactor = 0;
end
if phasres > 0
    nradiusphas = int(2 * nuesampling / phasres);
    if nradiusphas > nradius
        nradiusphas = nradius;
        warning( 'Warning: Not enough pts for requested phase res!')
    end
else
    nradiusphas = nradius;
end

nmax = 2^(nextpow2(nss+nradius)+zerofillingfactor);
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
% ifgphasre(1:length(xx)) = ifgacfwd(xx).*...
%         F_apoifg(apokind,nburstfwd,nradiusphas,xx);
    
    xx = nburstfwd:-1:nburstfwd-nradius+1;
    xx = xx(:);
    ifgre(end:-1:end-length(xx)+1) = 0.5*ifgacfwd(xx).*...
        F_apoifg(apokind,nburstfwd,nradius,xx);
    ifgphasre(end:-1:end-length(xx)+1) = ifgacfwd(xx).*...
        F_apophas(nburstfwd,nradiusphas,xx);
% ifgphasre(end:-1:end-length(xx)+1) = ifgacfwd(xx).*...
%         F_apoifg(apokind,nburstfwd,nradiusphas,xx);
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
specre = specre(1:nmax/2 + 1);
specim = specim(1:nmax/2 + 1);

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

