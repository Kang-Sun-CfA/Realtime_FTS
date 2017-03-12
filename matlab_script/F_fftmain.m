function [specrepc, phasre, phasim,nmax] = ...
    F_fftmain(ifgacfwd,nuesampling,nburstfwd,apokind,ssords,nss,nradius,phasres)
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
    disp([ 'OPDmax single-sided IFG:', num2str(nss / (2 * nuesampling))])
    if nradius > nss - 1
        nradius = nss - 1;
    end
    disp([ 'OPD available for phase calc:', num2str(nradius / (2 * nuesampling))])
    disp(['OPD requested for phase calc:', num2str(nradiusphas / (2 * nuesampling))])
    
    nmax = nmax * 2; % Apply additional zerofill of factor 2 for SS-IFG!!
    disp( ['FFT size: ', num2str(nmax)])
    deltanue = nuesampling / single(nmax / 2);
elseif ssords == 1
    % double-sided
    disp([ '   OPDmax: ',num2str((nradius / (2 * nuesampling)))])
    disp([ '   OPD requested for phase calc: ',...
        num2str((nradiusphas / (2 * nuesampling)))])
    disp( ['   FFT size: ', num2str(nmax)])
    deltanue = nuesampling / single(nmax / 2);
end
disp('   initializing FFT arrays ...')
ifgre = zeros(nmax,1,'single');
%     ifgim = zeros(nmax,1,'single');
ifgphasre = zeros(nmax,1,'single');
%     ifgphasim = zeros(nmax,1,'single');

if ssords == 0
    disp( '   before FFT, preparing SS-IFG')
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
    disp( '   before FFT, preparing DS-IFG')
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
