% This is the tutorial to use the working horse function,
% F_calmat(filename,config), where filename is the OPUS file, config is a
% structure containing input options.

% Written by Kang Sun on 2017/03/12

% system separator, "/" for mac and linux, "\" for windows
sfs = filesep;
% which EM27, we have ha and hb at harvard
observation = 'hb';
%% load OPUS
close all
for fileN = 1138;
    % example OPUS file located at Realtime_FTS/spectra/yymmdd/
    filename = ['..',sfs,'spectra',sfs,'160624',sfs,...
        observation,'20160624s0e00a.',sprintf('%0*d',4,fileN)];
    
    % initialize the input
    config = [];
    
    % which apodization.
    %        < 0: use the OPUS spectrum, no fft on the ifg
    %          1: boxcar
    %          2: triag
    %          3: Happ-Genzel
    %        4/5: Blackmann-Harris 3-term/4-term
    %      6/7/8: Norton-Beer weak/medium/strong
    config.apokind = 7;
    
    % lower and upper limit of output spectrum, in cm-1
    config.lowerlim = 5000;
    config.upperlim = 10000;
    
    % display processing message or not. I'd like it to be off
    config.display = false;
    
    % call the F_calmat function, you may time it by tic and toc:
    tic;calmat = F_calmat(filename,config);toc;
    
    % plot the results
    figure('color','w','unit','inch','position',[0 1 12 10])
    % ifg
    subplot(4,2,1:2)
    ax = plotyy(1:length(calmat.ifgC),calmat.ifgC,...
        1:length(calmat.ifgC),calmat.ifg_smooth);
    set(ax,'xlim',[1 length(calmat.ifgC)])
    legend('DC corrected ifg','DC of ifg');
    title(['ifgQuality = ',num2str(calmat.ifgQuality)])
    
    w2 = calmat.x;
    s2 = double(calmat.spec);
    % spectrum
    subplot(4,2,3:4)
    plot(w2,s2)
    subplot(4,2,5)
    plot(w2,s2);xlim([7765 8005])
    subplot(4,2,6)
    plot(w2,s2);xlim([6180 6260])
    subplot(4,2,7)
    plot(w2,s2);xlim([5880 5996])
    subplot(4,2,8)
    plot(w2,s2);xlim([6007 6145])
end
%%
[ifg,~,paramas] = ImportOpus(filename,'SampleInterferogram');
close all
figure('unit','inch','color','w','position',[-15 0 11 3.5])
subplot(1,2,1)
plot(1:length(calmat.ifgC),ifg,'k')
tmp = 0.003;
ylim([mean(calmat.ifg_smooth)-tmp,mean(calmat.ifg_smooth)+tmp])
set(gca,'linewidth',1,'xlim',[1,length(calmat.ifgC)])
xlabel('Number of samples')
ylabel('Signal intensity')
box off
title('Interferogram before correction','fontsize',14)
subplot(1,2,2)
plot(1:length(calmat.ifgC),calmat.ifgC,'k')
ylim([-tmp,+tmp])
set(gca,'linewidth',1,'xlim',[1,length(calmat.ifgC)],'ycolor','k')
box off
xlabel('Number of samples')
title('Interferogram after correction','fontsize',14)

%% plot ILS
% The function F_calmat_ils is essentially the same as F_calmat. The idea
% is that if you feed F_calmat with an ifg, it will give you a spectrum; if you
% feed F_calmat with a cosine, it will give an ILS. Zero filling
% of the cosine will oversample the ILS. It is also possible to add
% modulation loss and phase error this way, which could be more intuitive.

% extent of plotting the ILS, in wavenumber
ilsextent = 10; 
input_ils = [];

% need the metadata from just one OPUS file, any file should work
fileN = 1138;
filename = ['..',sfs,'spectra',sfs,'160624',sfs,...
        observation,'20160624s0e00a.',sprintf('%0*d',4,fileN)];
input_ils.filename = filename;

% wavenumber of the simulated cosine. The value should not matter, but
% fractions of nuesampling will make the ILS look better
input_ils.sim_wn = 1.579792480468750e+04/2;
% phase error. I'm not sure if it is correct though...
input_ils.phase_error = 1*5.412e-1;
% higher zero filling factor means more oversampling and more smooth ILS
input_ils.zerofillingfactor = 6;

apoarray = [1 6 7 8];
aponame = {'boxcar','NB weak','NB medium','NB strong'};
figure('color','w','unit','inch','position',[0 1 12 10])
hold on
for i = 1:length(apoarray)
input_ils.apokind = apoarray(i);
output_ils = F_calmat_ils(input_ils);
nmax = output_ils.nmax;
ilsx0 = output_ils.ilsx-input_ils.sim_wn;
int = ilsx0 >= -ilsextent & ilsx0 <= ilsextent;
ilsx0 = ilsx0(int);
ilsy0 = double(output_ils.ilsy(int)/sqrt(nmax));
plot(ilsx0,ilsy0)
end
legend(aponame);