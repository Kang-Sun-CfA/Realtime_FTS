clc;clear
cd ~/FTS/Realtime_FTS/
% Dir to save .mat intermediate data files
datadir = '/data/tempo1/Shared/kangsun/FTS_data/';
% HITRAN molecule parameter file
molparampath = '/home/kangsun/Courses/molparam.txt';

% fid = fopen([datadir,'av0.table']);
% C = cell2mat(textscan(fid,'%f%f%f%f','headerlines',1,'delimiter',' ',...
%     'multipledelimsasone',1));
% plot(C(:,1),C(:,2:4))
% %%
% vlow = C(:,1);
% [xgrid,ygrid] = meshgrid(vlow,[298 273 253]);
% h = surf(xgrid,ygrid,C(:,2:4)');set(h,'edgecolor','none')
% %%
% clc
% temp = 280;
% acia_low = interp2(xgrid,ygrid,C(:,2:4)',vlow,temp);
% plot(C(:,1),C(:,2:4),vlow,acia_low,'k')
%%%%%% Define resolution of common grid, in cm-1. Optical depth and solar
%%%%%% spectrum will be interpolated to this common grid
common_grid_resolution = 0.01;

common_grid = 7700:common_grid_resolution:8060;
vStart = min(common_grid);
vEnd = max(common_grid);

% load a high profile
fid = fopen('/data/tempo2/ggonzale/GEOCAPE-TOOL/geocape_data/TES_Data/AtmProfiles_Data_Run551C_Seq0028_Scan002.asc');
C = cell2mat(textscan(fid,'%f%f%f%f%*[^\n]','headerlines',10,'multipledelimsasone',1));
fclose(fid);

load([datadir,'pt0624.mat'],'T_prof','H_prof','H2O_prof','P_prof')

% plot(T_prof,H_prof,C(:,4),C(:,3)*1e3)
pTES = C(:,2);tTES = C(:,4);hTES = C(:,3)*1000; h2oTES = C(:,1)*0;

int = pTES < min(P_prof);
P_prof = [flipud(pTES(int));P_prof];
H_prof = [flipud(hTES(int));H_prof];
T_prof = [flipud(tTES(int));T_prof];
H2O_prof = [flipud(h2oTES(int));H2O_prof];

plot(T_prof,H_prof)


nlevel = length(P_prof);
nlayer = nlevel-1;
input_LBL = [];
input_LBL.P_level = P_prof;
input_LBL.P_layer = P_prof(1:end-1)+diff(P_prof)/2;
input_LBL.T_layer = interp1(P_prof,T_prof,input_LBL.P_layer);
input_LBL.Z_level = H_prof;
input_LBL.Z_layer = interp1(P_prof,H_prof,input_LBL.P_layer);
input_LBL.taumol = 'O2';
input_LBL.profile = ones(nlayer,1)*0.2095;
input_LBL.common_grid = common_grid;
input_LBL.molparampath = molparampath;
input_LBL.datadir = datadir;
input_LBL.which_CIA = {'gfit','sao','mate'};% keep lower case
output_LBL = F_line_by_line(input_LBL);

input_LBL.taumol = 'H2O';
input_LBL.profile = interp1(P_prof,H2O_prof,input_LBL.P_layer)/0.622;
output_LBL_H2O = F_line_by_line(input_LBL);

%% load OPUS
clc
wStart = 7765;
wEnd = 8005;
apokind_all = 7;
observation = 'hb';
fileN = 1856;
filename = ['/data/tempo1/Shared/kangsun/FTS_data/EM27_data/em27spectra/160624/',observation,'20160624s0e00a.',sprintf('%0*d',4,fileN)];
config = [];
config.apokind = apokind_all;
% config.refraction = 1;
% config.pathpt = 'C:\\FTS\\em27spectra\\160624\\pt';
config.lowerlim = wStart;
config.upperlim = wEnd;
config.localt = 0.22;
config.loclat = 39.766667;
config.loclon = -86.15;
calmat = F_calmat(filename,config);
w2 = calmat.x;
s2 = double(calmat.spec);
plot(w2,s2)
%% simulate ILS
ilsextent = 10; % wavenumber
input_ils = [];
input_ils.filename = filename;
input_ils.sim_wn = 1.579792480468750e+04/2;
input_ils.phase_error = 0*5.412e-1;
input_ils.apokind = apokind_all;
input_ils.zerofillingfactor = 6;
output_ils = F_calmat_ils(input_ils);
%
nmax = output_ils.nmax;
ilsx0 = output_ils.ilsx-input_ils.sim_wn;
int = ilsx0 >= -ilsextent & ilsx0 <= ilsextent;
ilsx0 = ilsx0(int);
ilsy0 = double(output_ils.ilsy(int)/sqrt(nmax));
ilsx = -ilsextent:common_grid_resolution:ilsextent;
ilsy = interp1(ilsx0,ilsy0,ilsx,'linear','extrap');
% ilsy(isnan(ilsy)) = 0;
ilsy = ilsy/abs(sum(ilsy));
plotyy(ilsx0,ilsy0,ilsx,ilsy)
%% load solar data and convolve/resample to common grd
load([datadir, 'FTS_solar.mat'],'llwaven','lltranswaven')
int = llwaven >= vStart & llwaven <= vEnd;
llwave = llwaven(int);lltrans = lltranswaven(int);
% Bring solar spectrum to common_grid
input_nlinfit = [];
input_nlinfit.ss = ...
    F_conv_interp(llwave,lltrans,...
    2.5*common_grid_resolution,common_grid);
plot(llwave,lltrans,common_grid,input_nlinfit.ss)
%% prepare fitting database of absorptions
input_nlinfit.s1 = nan(3,length(common_grid));
input_nlinfit.s1(1,:) = sum(output_LBL.dtau_cntm_sao,2);% CIA
input_nlinfit.s1(2,:) = squeeze(sum(output_LBL.dtau_layer(1,:,:),3));% O2
input_nlinfit.s1(3,:) = squeeze(sum(output_LBL_H2O.dtau_layer(1,:,:),3));% H2O

input_nlinfit.ils = ilsy;
input_nlinfit.w2 = w2;
input_nlinfit.w1 = common_grid;
%% prepare optional inputs into the nlinfit
clc
% remove all leftover optional fields in input_nlinfit
optional_fields = {'shift_gas','shift_sun','shift_together','scaling',...
    'tilt','zlo'};
for ifield = 1:length(optional_fields)
    if isfield(input_nlinfit,optional_fields{ifield})
        input_nlinfit = rmfield(input_nlinfit,optional_fields{ifield});
    end
end
input_nlinfit.shift_gas = 0;
% input_nlinfit.shift_sun = 0;
input_nlinfit.shift_together = 0;
% input_nlinfit.scaling = 1;
% input_nlinfit.tilt = 0;
% input_nlinfit.zlo = 0;

coeff0 = [0 .2 0 0 16 10 0];

input_nlinfit.s1(1,:) = sum(output_LBL.dtau_cntm_sao,2);% CIA
tic
[coeff_sao, R_sao] = nlinfit(input_nlinfit,s2,@F_forward,coeff0);
toc

input_nlinfit.s1(1,:) = sum(output_LBL.dtau_cntm_gfit,2);% CIA
tic
[coeff_gfit, R_gfit] = nlinfit(input_nlinfit,s2,@F_forward,coeff0);
toc

input_nlinfit.s1(1,:) = sum(output_LBL.dtau_cntm_mate,2);% CIA
tic
[coeff_mate, R_mate] = nlinfit(input_nlinfit,s2,@F_forward,coeff0);
toc
%%
close all
figure('unit','inch','position',[0 1 12 5])
subplot(2,1,1)
input_nlinfit.s1(1,:) = sum(output_LBL.dtau_cntm_sao,2);% CIA
hold on
plot(w2,s2,'k',w2,F_forward(coeff_sao,input_nlinfit))
input_nlinfit.s1(1,:) = sum(output_LBL.dtau_cntm_gfit,2);% CIA
plot(w2,F_forward(coeff_gfit,input_nlinfit))
input_nlinfit.s1(1,:) = sum(output_LBL.dtau_cntm_mate,2);% CIA
plot(w2,F_forward(coeff_mate,input_nlinfit))
legend('FTS','SAO','GFIT','MATE')
ylabel('Spectra')
subplot(2,1,2)
plot(w2,R_sao,w2,R_gfit,w2,R_mate)
legend('SAO','GFIT','MATE')
xlabel('Wavenumber')
ylabel('Residuals')
%%
export_fig 'fitting_residuals_ani_127.pdf'
%%
close all
figure('unit','inch','position',[0 1 10 5])
plot(common_grid,sum(output_LBL.dtau_cntm_sao,2),...
    common_grid,sum(output_LBL.dtau_cntm_gfit,2),...
    common_grid,sum(output_LBL.dtau_cntm_mate,2))
xlabel('Wavenumber')
ylabel('Optical depth')
legend('SAO','GFIT','MATE')
%%
export_fig 'CIA_comparison_ani_127.pdf'
%%
plot(w2,F_forward(coeff0,input_nlinfit),w2,s2,'k',w2,F_forward(coeff,input_nlinfit))
plot(w2,F_forward(coeff0,input_nlinfit))