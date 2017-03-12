% functionpath = 'D:\Research_CfA\FTS\Matlab_scripts';
% OPD = 1.8;
% FOV = 2.36e-3; % in rad
% phase_error = 0*5.412e-3;
% modulation_loss = 1-9.871e-1;
% apd = 2;

clc
% close all
cd(functionpath)
S_input_parameters_LV;

C = [1 0 0 0 0;
    0.384093 -0.087577 0.703484 0 0;
    0.152442 -0.136176 0.983734 0 0;
    0.045335 0         0.554883 0 0.399782];

ilsx = -ilsextent:common_grid_resolution:ilsextent;
ILS_input = nan(length(window_list),length(ilsx));
for iwin = 1:length(window_list)
    
    target_gas = window_list(iwin).target_gas;
    windowID = window_list(iwin).windowID;
    [~, ~, wStart, wEnd,~,~] = ...
        F_define_windows(target_gas,windowID);
    Fl = mean([wStart wEnd]);
    
    Fs = Fl*8;                   % Sampling frequency
    T = 1/Fs;                     % Sampling period
    L = floor(OPD*Fs);                     % Length of signal
    n = 2^nextpow2(L)*32;         % length of time domain signal, padding 0
    t = (0:L-1)*T;                % Time vector
    U = linspace(0,1,length(t));
    FN = zeros(size(t));
    for i = 0:4;
        FN = FN+C(apd+1,i+1)*(1-U.^2).^i;
    end
    
    aperture_x = pi*FOV^2*sinc(Fl*t*pi*FOV^2/2/pi);
    ifg = cos(2*pi*Fl*t-phase_error)...
        .*linspace(1,1-modulation_loss,length(t))...
        .*FN.*aperture_x;
    % ifg = aperture_x;
    
    Y = fft(ifg,n);
    P2 = real(Y/n);
    P1 = P2(:,1:n/2+1);
    P1(:,2:end-1) = 2*P1(:,2:end-1);
    xdata = 0:(Fs/n):(Fs/2-Fs/n);
    ydata1 = P1(1,1:n/2)/max(abs(P1(1,1:n/2)));
    ILS_input(iwin,:) = interp1(xdata-Fl,ydata1,ilsx);
%     figure
%     plot(xdata,ydata1,'.')
%     xlim([Fl-ilsextent,Fl+ilsextent])
end