% overhaul the fit CIA code to use absco tables
% updated by Kang Sun on 2017/10/23

clear;clc;close all
% system separator, "/" for mac and linux, "\" for windows
sfs = filesep;

% Load fitting window information and pre-calculated optical depth
profile_date = '20160624';% yyyymmdd, 20170000 for US standard atmosphere
% profile_date = '20170217';% yyyymmdd, 20170000 for US standard atmosphere

profile_hour = '0000';% hhmm, 0000 for US standard atmosphere
window_list_fn = ['..',sfs,'spectroscopy',sfs,'window_list_',...
    profile_date,'_',profile_hour,'_cia.mat'];
inp_tau = [];

%!!!!!!!!!!!!!!!! local dir saving ABSCO tables !!!!!!!!!!!!!!!!!!!!!!
inp_tau.absco_dir = '/data/tempo1/Shared/kangsun/FTS_data/HITRAN/';
%!!!!!!!!!!!!!!!! YOU HAVE TO SPECIFY THIS !!!!!!!!!!!!!!!!!!!!!!!!!!!

inp_tau.CIA_dir = ['..',sfs,'spectroscopy',sfs];
inp_tau.which_CIA = {'sao','gfit','mate','hitran','koenis'};

inp_tau.profile_fn = ['..',sfs,'profiles',sfs,'hb',profile_date,...
    '.map'];
inp_tau.lowest_possible_Psurf = 980;
inp_tau.nsublayer = 1;
inp_tau.common_grid_resolution = 0.01;
inp_tau.surface_layer_P = 990;

% end of input

window_list = F_absco2tau_flexible(inp_tau);
save(window_list_fn,'window_list')
%%
load(window_list_fn)
% %% add HDO and potentially other non-absco molecules
% % which window to add more molecules?
% iwin = 2;
% % load US standard atmosphere
% fid = fopen('../profiles/profiles_20170000_0000.dat');
% C = cell2mat(textscan(fid,'%f%f%f%f%f%f%f','delimiter',' ','multipledelimsasone',1,'headerlines',11));
% 
% % Dir to save .mat intermediate data files
% datadir = '../spectroscopy/';
% % HITRAN molecule parameter file
% molparampath = '/home/kangsun/Courses/molparam.txt';
% 
% P_prof = C(:,1);
% T_prof = C(:,3);
% H_prof = C(:,2);
% H2O_prof = C(:,4);
% 
% nlevel = length(P_prof);
% nlayer = nlevel-1;
% input_LBL = [];
% input_LBL.P_level = P_prof;
% input_LBL.P_layer = P_prof(1:end-1)+diff(P_prof)/2;
% input_LBL.T_layer = interp1(P_prof,T_prof,input_LBL.P_layer);
% input_LBL.Z_level = H_prof;
% input_LBL.Z_layer = interp1(P_prof,H_prof,input_LBL.P_layer);
% input_LBL.common_grid = window_list(iwin).common_grid;
% input_LBL.taumol = 'H2O';
% input_LBL.profile = interp1(P_prof,H2O_prof,input_LBL.P_layer);
% input_LBL.molparampath = molparampath;
% input_LBL.datadir = datadir;
% input_LBL.which_CIA = {'gfit','sao','mate'};% keep lower case
% 
% output_LBL_H2O = F_line_by_line(input_LBL);
% %% compare optical depth of H2O from absco and lbl
% % I found that the 1.06 um window does not have HDO.
% close
% subplot(2,1,1)
% plot(window_list(iwin).common_grid,squeeze(sum(output_LBL_H2O.dtau_layer(1,:,:),3)))
% subplot(2,1,2)
% plot(window_list(iwin).common_grid,squeeze(sum(output_LBL_H2O.dtau_layer(2,:,:),3)))
% %% check hdo
% load('../spectroscopy/FTS_5000to12000.mat');
% int = lines.moleculeNumber == 1 & lines.isotopologueNumber == 4;
% plot(lines.transitionWavenumber(int),lines.lineIntensity(int),'.')
% xlim([7700 8050])
% xlim([9100 9700])
% %% test the effect of broadening. Results: it does show difference in tau, but not so much in residual. Disappointing
% clear;clc;close all
% % system separator, "/" for mac and linux, "\" for windows
% sfs = filesep;
% 
% % Load fitting window information and pre-calculated optical depth
% profile_date = '20160624';% yyyymmdd, 20170000 for US standard atmosphere
% % profile_date = '20170217';% yyyymmdd, 20170000 for US standard atmosphere
% 
% profile_hour = '0000';% hhmm, 0000 for US standard atmosphere
% window_list_fn = ['..',sfs,'spectroscopy',sfs,'window_list_',...
%     profile_date,'_',profile_hour,'_cia.mat'];
% inp_tau = [];
% 
% %!!!!!!!!!!!!!!!! local dir saving ABSCO tables !!!!!!!!!!!!!!!!!!!!!!
% inp_tau.absco_dir = '/data/tempo1/Shared/kangsun/FTS_data/HITRAN/';
% %!!!!!!!!!!!!!!!! YOU HAVE TO SPECIFY THIS !!!!!!!!!!!!!!!!!!!!!!!!!!!
% 
% inp_tau.CIA_dir = ['..',sfs,'spectroscopy',sfs];
% inp_tau.which_CIA = {'sao','gfit','mate','hitran','koenis'};
% 
% inp_tau.profile_fn = ['..',sfs,'profiles',sfs,'hb',profile_date,...
%     '.map'];
% inp_tau.lowest_possible_Psurf = 980;
% inp_tau.nsublayer = 1;
% inp_tau.common_grid_resolution = 0.01;
% inp_tau.surface_layer_P = 990;
% 
% % end of input
% 
% window_list1 = F_absco2tau_flexible(inp_tau);
% inp_tau.if_broadening = true;
% window_list2 = F_absco2tau_flexible(inp_tau);
% %%
% iwin = 2;
% plot(window_list1(iwin).common_grid,window_list1(iwin).tau_struct.H2O.Tau_sum,...
%     window_list2(iwin).common_grid,window_list2(iwin).tau_struct.H2O.Tau_sum)
% %%
% plot(window_list1(1).common_grid,window_list1(1).tau_struct.O2.Tau_sum,...
%     window_list2(1).common_grid,window_list2(1).tau_struct.O2.Tau_sum)
%% 
common_grid_resolution = mean(diff(window_list(1).common_grid));
% Half length of the ILS, the longer the more accurate and slower
ilsextent = 10; % wavenumber
% find an OPUS file with the same ILS as the measurements
ils_fn = '/data/tempo1/Shared/kangsun/FTS_data/EM27_data/em27spectra/160624/hb20160625s0e00a.0019';

input_ils = [];
input_ils.filename = ils_fn;
input_ils.sim_wn = 1.579792480468750e+04/2;
input_ils.phase_error = 0*5.412e-1;
input_ils.apokind = 7; % should be consistent with apokind in config
input_ils.zerofillingfactor = 6;
output_ils = F_calmat_ils(input_ils);
nmax = output_ils.nmax;
ilsx0 = output_ils.ilsx-input_ils.sim_wn;
int = ilsx0 >= -ilsextent & ilsx0 <= ilsextent;
ilsx0 = ilsx0(int);
ilsy0 = double(output_ils.ilsy(int)/sqrt(nmax));
ilsx = -ilsextent:common_grid_resolution:ilsextent;
ilsy = interp1(ilsx0,ilsy0,ilsx,'linear','extrap');
ilsy = ilsy/abs(sum(ilsy));

% initialize the input
config = [];

config.apokind = 7;

% lower and upper limit of output spectrum, in cm-1
config.lowerlim = 5000;
config.upperlim = 10000;

% display processing message or not. I'd like it to be off
config.display = false;

config.ifamf = true;

% coordinate and altitude (km) of measurement location
config.loclat = 42.377162;
config.loclon = -71.113461;
config.localt = 0.22;

% time difference compared to UTC in hour
config.timezone = -5;
filename = '/data/tempo1/Shared/kangsun/FTS_data/EM27_data/em27spectra/160624/hb20160625s0e00a.0019';

calmat = F_calmat(filename,config);
%%
common_grid_resolution = mean(diff(window_list(1).common_grid));
% Half length of the ILS, the longer the more accurate and slower
ilsextent = 10; % wavenumber
% find an OPUS file with the same ILS as the measurements
ils_fn = ['..',sfs,'spectra',sfs,'160624',sfs,'hb20160624s0e00a.1138'];
input_ils = [];
input_ils.filename = ils_fn;
input_ils.sim_wn = 1.579792480468750e+04/2;
input_ils.phase_error = 0*5.412e-1;
input_ils.apokind = 7; % should be consistent with apokind in config
input_ils.zerofillingfactor = 6;
output_ils = F_calmat_ils(input_ils);
nmax = output_ils.nmax;
ilsx0 = output_ils.ilsx-input_ils.sim_wn;
int = ilsx0 >= -ilsextent & ilsx0 <= ilsextent;
ilsx0 = ilsx0(int);
ilsy0 = double(output_ils.ilsy(int)/sqrt(nmax));
ilsx = -ilsextent:common_grid_resolution:ilsextent;
ilsy = interp1(ilsx0,ilsy0,ilsx,'linear','extrap');
ilsy = ilsy/abs(sum(ilsy));

% initialize the input
config = [];

config.apokind = 7;

% lower and upper limit of output spectrum, in cm-1
% config.lowerlim = 5000;
% config.upperlim = 10000;

% display processing message or not. I'd like it to be off
config.display = false;

config.ifamf = true;

% coordinate and altitude (km) of measurement location
config.loclat = 39.767;
config.loclon = -86.15;
config.localt = 0.22;

% time difference compared to UTC in hour
config.timezone = -5;
filename = '/data/tempo1/Shared/kangsun/FTS_data/EM27_data/em27spectra/160624/hb20160625s0e00a.0019';
calmat = F_calmat(filename,config);
%%
Results = cell(length(window_list),2);
load(['..',sfs,'spectroscopy',sfs,'FTS_solar.mat'],'llwaven','lltranswaven')
%%
clc
P_surface = 1020;
clc; close all
input_nlinfit = [];
input_nlinfit.ils = ilsy;
optional_fields = {'shift_gas','shift_sun','shift_together','scaling','tilt','zlo'};
optional_fields_ap = [0,           0,          0,               1,       0,     0];
% the following index combined should cover
% 1:length(optional_fields)
input_fields_index = [ 3];% index as input
fit_fields_index = [1 2 4 5 6];% index included in fitting
input_fields = optional_fields(input_fields_index);
input_value = optional_fields_ap(input_fields_index);
for i = 1:length(input_value)
    input_nlinfit.(input_fields{i}) = input_value(i);
end
if ~exist('ss','var')
    ss = cell(length(window_list),1);
end

CC = lines(5);
window_list(2).wRange = [7765, 8035];
window_list(1).wRange = [9200, 9590];

addpath('~/matlab functions/export_fig/')
which_CIA = 'GFIT';

iwin = 1;%length(window_list)
target_gas = window_list(iwin).target_gas;
windowID = window_list(iwin).windowID;
common_grid = window_list(iwin).common_grid;

vStart = window_list(iwin).vRange(1);
vEnd = window_list(iwin).vRange(2);

wStart = window_list(iwin).wRange(1);
wEnd = window_list(iwin).wRange(2);

interval = calmat.x >= wStart & calmat.x <= wEnd;

w2 = calmat.x(interval);
s2 = double(calmat.spec(interval));
w2 = w2(:);s2 = s2(:)/max(s2);

if isempty(ss{iwin})
    int = llwaven >= vStart & llwaven <= vEnd;
    llwave = llwaven(int);lltrans = lltranswaven(int);
    ss{iwin} = F_conv_interp(llwave,lltrans,...
        2.5*common_grid_resolution,common_grid);
end
input_nlinfit.ss = ss{iwin};

% this is very similar to mol_for_fit, except O2 CIA
fit_spec_names = fieldnames(window_list(iwin).tau_struct);
n_fit_spec = length(fit_spec_names);
input_nlinfit.s1 = nan(n_fit_spec,length(common_grid));
for i_fit_spec = 1:n_fit_spec
    tmp_tau_prof = window_list(iwin).tau_struct.(fit_spec_names{i_fit_spec}).Tau_sum;
    tmp_tau_sfc = window_list(iwin).tau_struct.(fit_spec_names{i_fit_spec}).surface_layer_sigma/1e4... % absorption cross section in m2
        *(P_surface-window_list(iwin).tau_struct.(fit_spec_names{i_fit_spec}).surface_layer_top_P)*100 ... % surface layer thickness in Pa
        /(0.029/6.02e23*9.8) ... % air molecular weight and g
        *window_list(iwin).tau_struct.(fit_spec_names{i_fit_spec}).surface_layer_VMR;
    input_nlinfit.s1(i_fit_spec,:) = double(tmp_tau_prof+tmp_tau_sfc);
end

n_fit_spec = n_fit_spec+1;
tau_CIA = double(window_list(iwin).tau_struct.O2.(['CIA_',which_CIA]));
%     tau_CIA = double(window_list(iwin).tau_struct.O2.CIA_GFIT);
input_nlinfit.s1 = cat(1,input_nlinfit.s1,tau_CIA);

input_nlinfit.ils = ilsy;
input_nlinfit.w2 = w2;
input_nlinfit.w1 = common_grid;

%             input_nlinfit.shift_gas = 0;
%             % input_nlinfit.shift_sun = 0;
input_nlinfit.shift_together = 0;
%             % input_nlinfit.scaling = 1;
%             % input_nlinfit.tilt = 0;
%             % input_nlinfit.zlo = 0;

if iwin == 1
    input_nlinfit.coeff_bound = [-inf, -inf, -inf, -inf, -inf, -inf, -inf, 16, -inf;
        inf, inf, inf, inf, inf, inf, inf, 18, inf];
end
coeff0 = [optional_fields_ap(fit_fields_index) ones(1,n_fit_spec)*16.878];
[coeff, R] = nlinfit(input_nlinfit,s2,@F_forward,coeff0);

f_gfit = F_forward(coeff,input_nlinfit);
coeff_gfit = coeff;
R_gfit = R;

which_CIA = 'HITRAN';
input_nlinfit.s1(end,:) = double(window_list(iwin).tau_struct.O2.(['CIA_',which_CIA]));
if iwin == 1
    input_nlinfit.coeff_bound = [-inf, -inf, -inf, -inf, -inf, -inf, -inf, 16, -inf;
        inf, inf, inf, inf, inf, inf, inf, 18, inf];
end
coeff0 = [optional_fields_ap(fit_fields_index) ones(1,n_fit_spec)*16.878];
[coeff, R] = nlinfit(input_nlinfit,s2,@F_forward,coeff0);

f_theory = F_forward(coeff,input_nlinfit);
coeff_theory = coeff;
R_theory = R;

which_CIA = 'Koenis';
input_nlinfit.s1(end,:) = double(window_list(iwin).tau_struct.O2.(['CIA_',which_CIA]));
if iwin == 1
    input_nlinfit.coeff_bound = [-inf, -inf, -inf, -inf, -inf, -inf, -inf, 16, -inf;
        inf, inf, inf, inf, inf, inf, inf, 18, inf];
end
coeff0 = [optional_fields_ap(fit_fields_index) ones(1,n_fit_spec)*16.878];
[coeff, R] = nlinfit(input_nlinfit,s2,@F_forward,coeff0);

f_koenis = F_forward(coeff,input_nlinfit);
coeff_koenis = coeff;
R_koenis = R;

%%
stretchx = 1.1;
stretchy = 1.15;
close all
figure('unit','inch','color','w','position',[0 0 12 7])
subplot(3,2,2)
plot(w2,s2,'k')
ylim([0 1.05])
xlim([9200 9600])
set(gca,'linewidth',1,'box','off')
set(gca,'xcolor','none')
pos = get(gca,'position');
set(gca,'position',[pos(1:2) pos(3)*stretchx pos(4)*stretchy],'xgrid','on')
title('(d) $X^3\Sigma_g^-(v=0) \rightarrow a^1\Delta_g(v=1)$: observation','interpreter','latex')

subplot(3,2,4)
hold on
input_nlinfit.coeff_bound = [-inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf;
        inf, inf, inf, inf, inf, inf, inf, inf, inf];
plot(w2,F_forward([coeff_gfit(1:5) 0 coeff_gfit(7) 0 0],input_nlinfit))
plot(w2,F_forward([coeff_gfit(1:5) 0 0 0 coeff_gfit(9)],input_nlinfit))
plot(w2,F_forward([coeff_gfit(1:5) coeff_gfit(6) 0 0 0],input_nlinfit))
plot(w2,F_forward([coeff_gfit(1:5) 0 0 coeff_gfit(8) 0],input_nlinfit))
plot(w2,F_forward([coeff_gfit(1:5) 0 0 0 0],input_nlinfit),'k')
ylim([0 1.05])
xlim([9200 9600])
set(gca,'linewidth',1,'box','off')
set(gca,'xcolor','none')
pos = get(gca,'position');
set(gca,'position',[pos(1:2) pos(3)*stretchx pos(4)*stretchy],'xgrid','on')
hleg = legend('H_2O','O_2 CIA','O_2','CO_2','Solar');
set(hleg,'location','southeast','box','off')
title('(e) $X^3\Sigma_g^-(v=0) \rightarrow a^1\Delta_g(v=1)$: absorbers','interpreter','latex')

subplot(3,2,6)
h = plot(w2,R_gfit,'k',w2,R_theory,w2,R_koenis);%,w2,s2-f_gfit)
ylim([-0.05 0.05])
xlim([9200 9600])
set(gca,'linewidth',1,'box','off')
pos = get(gca,'position');
set(gca,'position',[pos(1:2) pos(3)*stretchx pos(4)*stretchy],'xgrid','on')
hleg = legend('GGG2014','Theory','Experiment');
set(hleg,'location','south','orientation','hori','box','off')

xlabel('Frequency [cm$^{-1}$]','interpreter','latex')
title('(f) $X^3\Sigma_g^-(v=0) \rightarrow a^1\Delta_g(v=1)$: fitting residuals','interpreter','latex')

%% 1.27 um
P_surface = 1020;
clc; 
input_nlinfit = [];
input_nlinfit.ils = ilsy;
optional_fields = {'shift_gas','shift_sun','shift_together','scaling','tilt','zlo'};
optional_fields_ap = [0,           0,          0,               1,       0,     0];
% the following index combined should cover
% 1:length(optional_fields)
input_fields_index = [ 3];% index as input
fit_fields_index = [1 2 4 5 6];% index included in fitting
input_fields = optional_fields(input_fields_index);
input_value = optional_fields_ap(input_fields_index);
for i = 1:length(input_value)
    input_nlinfit.(input_fields{i}) = input_value(i);
end
if ~exist('ss','var')
    ss = cell(length(window_list),1);
end

CC = lines(5);
window_list(2).wRange = [7765, 8035];
window_list(1).wRange = [9200, 9590];

iwin = 2;%length(window_list)
target_gas = window_list(iwin).target_gas;
windowID = window_list(iwin).windowID;
common_grid = window_list(iwin).common_grid;

vStart = window_list(iwin).vRange(1);
vEnd = window_list(iwin).vRange(2);

wStart = window_list(iwin).wRange(1);
wEnd = window_list(iwin).wRange(2);

interval = calmat.x >= wStart & calmat.x <= wEnd;

w2 = calmat.x(interval);
s2 = double(calmat.spec(interval));
w2 = w2(:);s2 = s2(:)/max(s2);

if isempty(ss{iwin})
    int = llwaven >= vStart & llwaven <= vEnd;
    llwave = llwaven(int);lltrans = lltranswaven(int);
    ss{iwin} = F_conv_interp(llwave,lltrans,...
        2.5*common_grid_resolution,common_grid);
end
input_nlinfit.ss = ss{iwin};

% this is very similar to mol_for_fit, except O2 CIA
fit_spec_names = fieldnames(window_list(iwin).tau_struct);
n_fit_spec = length(fit_spec_names);
input_nlinfit.s1 = nan(n_fit_spec,length(common_grid));
for i_fit_spec = 1:n_fit_spec
    tmp_tau_prof = window_list(iwin).tau_struct.(fit_spec_names{i_fit_spec}).Tau_sum;
    tmp_tau_sfc = window_list(iwin).tau_struct.(fit_spec_names{i_fit_spec}).surface_layer_sigma/1e4... % absorption cross section in m2
        *(P_surface-window_list(iwin).tau_struct.(fit_spec_names{i_fit_spec}).surface_layer_top_P)*100 ... % surface layer thickness in Pa
        /(0.029/6.02e23*9.8) ... % air molecular weight and g
        *window_list(iwin).tau_struct.(fit_spec_names{i_fit_spec}).surface_layer_VMR;
    input_nlinfit.s1(i_fit_spec,:) = double(tmp_tau_prof+tmp_tau_sfc);
end

n_fit_spec = n_fit_spec+1;
% tau_CIA = double(window_list(iwin).tau_struct.O2.(['CIA_',which_CIA]));
tau_CIA = double(window_list(iwin).tau_struct.O2.CIA_GFIT);
input_nlinfit.s1 = cat(1,input_nlinfit.s1,tau_CIA);

input_nlinfit.ils = ilsy;
input_nlinfit.w2 = w2;
input_nlinfit.w1 = common_grid;

%             input_nlinfit.shift_gas = 0;
%             % input_nlinfit.shift_sun = 0;
input_nlinfit.shift_together = 0;
%             % input_nlinfit.scaling = 1;
%             % input_nlinfit.tilt = 0;
%             % input_nlinfit.zlo = 0;


coeff0 = [optional_fields_ap(fit_fields_index) ones(1,n_fit_spec)*16.878];
[coeff, R] = nlinfit(input_nlinfit,s2,@F_forward,coeff0);

f_gfit = F_forward(coeff,input_nlinfit);
coeff_gfit = coeff;
R_gfit = R;

which_CIA = 'HITRAN';
input_nlinfit.s1(end,:) = double(window_list(iwin).tau_struct.O2.(['CIA_',which_CIA]));
if iwin == 1
    input_nlinfit.coeff_bound = [-inf, -inf, -inf, -inf, -inf, -inf, -inf, 16, -inf;
        inf, inf, inf, inf, inf, inf, inf, 18, inf];
end
coeff0 = [optional_fields_ap(fit_fields_index) ones(1,n_fit_spec)*16.878];
[coeff, R] = nlinfit(input_nlinfit,s2,@F_forward,coeff0);

f_theory = F_forward(coeff,input_nlinfit);
coeff_theory = coeff;
R_theory = R;

which_CIA = 'Mate';
input_nlinfit.s1(end,:) = double(window_list(iwin).tau_struct.O2.(['CIA_',which_CIA]));
if iwin == 1
    input_nlinfit.coeff_bound = [-inf, -inf, -inf, -inf, -inf, -inf, -inf, 16, -inf;
        inf, inf, inf, inf, inf, inf, inf, 18, inf];
end
coeff0 = [optional_fields_ap(fit_fields_index) ones(1,n_fit_spec)*16.878];
[coeff, R] = nlinfit(input_nlinfit,s2,@F_forward,coeff0);

f_mate = F_forward(coeff,input_nlinfit);
coeff_mate = coeff;
R_mate = R;
%%
subplot(3,2,1)
plot(w2,s2/coeff_gfit(3),'k')
ylim([0 1.05])
xlim([7765, 8035])
set(gca,'linewidth',1,'box','off')
set(gca,'xcolor','none')
pos = get(gca,'position');
set(gca,'position',[pos(1:2) pos(3)*stretchx pos(4)*stretchy],'xgrid','on')
title('(a) $X^3\Sigma_g^-(v=0) \rightarrow a^1\Delta_g(v=0)$: observation','interpreter','latex')

subplot(3,2,3)
hold on
plot(w2,F_forward([coeff_gfit(1:5) 0 coeff_gfit(7) 0],input_nlinfit)/coeff_gfit(3))
plot(w2,F_forward([coeff_gfit(1:5) 0 0 coeff_gfit(8)],input_nlinfit)/coeff_gfit(3))
plot(w2,F_forward([coeff_gfit(1:5) coeff_gfit(6) 0 0],input_nlinfit)/coeff_gfit(3))
plot(w2,F_forward([coeff_gfit(1:5) 0 0 0],input_nlinfit)/coeff_gfit(3),'k')
ylim([0 1.05])
xlim([7765, 8035])
set(gca,'linewidth',1,'box','off')
hold off
set(gca,'xcolor','none')
pos = get(gca,'position');
set(gca,'position',[pos(1:2) pos(3)*stretchx pos(4)*stretchy],'xgrid','on')
hleg = legend('H_2O','O_2 CIA','O_2','Solar');
set(hleg,'location','southeast','box','off')
title('(b) $X^3\Sigma_g^-(v=0) \rightarrow a^1\Delta_g(v=0)$: absorbers','interpreter','latex')


subplot(3,2,5)
h=plot(w2,R_gfit/coeff_gfit(3),'k',w2,R_theory/coeff_gfit(3),w2,R_mate/coeff_gfit(3));%,w2,s2-f_gfit)
xlim([7765, 8035])
ylim([-0.05 0.05])
set(gca,'linewidth',1,'box','off')
pos = get(gca,'position');
set(gca,'position',[pos(1:2) pos(3)*stretchx pos(4)*stretchy],'xgrid','on')
hleg = legend('GGG2014','Theory','Experiment');
set(hleg,'location','south','orientation','hori','box','off')
% uistack(h(1),'top')
% uistack(h(2),'top')

xlabel('Frequency [cm$^{-1}$]','interpreter','latex')
title('(c) $X^3\Sigma_g^-(v=0) \rightarrow a^1\Delta_g(v=0)$: fitting residuals','interpreter','latex')
%%
export_fig('../figures/CIA_validation.pdf')
%%

   
%     close
    figure('unit','inch','color','w','position',[5 1 7 7])
    ax1 = subplot(2,1,1);
    bsln = polyval([coeff(4), coeff(3)],w2-mean(w2));
    h = plot(w2,s2./bsln-coeff(5),'k',w2,F_forward(coeff,input_nlinfit)./bsln-coeff(5),'linewidth',.5);
    set(h(2),'color',CC(3,:))
    hold on
%     tmp = conv(input_nlinfit.ss,ilsy/sum(ilsy),'same');
%     tmp = interp1(input_nlinfit.w1,tmp,w2,'linear','extrap');
%     plot(w2,tmp)
    for is1 = size(input_nlinfit.s1,1)
        tmp = conv(exp(-input_nlinfit.s1(is1,:)*coeff(5+is1)),ilsy/sum(ilsy),'same');
    tmp = interp1(input_nlinfit.w1,tmp,w2,'linear','extrap');
    plot(w2,tmp,'linewidth',1.5)
    end
    set(gca,'xticklabel',[])
    title(which_CIA)
    ylim([-0.05 1.05])
    hleg = legend('Observation','Fitting','CIA transmission');set(hleg,'box','off')
    if iwin == 1;set(hleg,'location','east');end
    
    ax2 = subplot(2,1,2);
    h = plot(w2,(s2-F_forward(coeff,input_nlinfit))./bsln,'k','linewidth',1.);
    pos = get(ax1,'position');
    set(ax1,'linewidth',1,'position',[pos(1) pos(2)-0.2 pos(3) pos(4)+0.2],...
        'xlim',[wStart wEnd])
    
    pos = get(ax2,'position');
    set(ax2,'linewidth',1,'position',[pos(1) pos(2) pos(3) pos(4)-0.15],...
        'xlim',[wStart wEnd])
    xlabel('Wavenumber')
    ylim([-0.05 0.05])
    hleg = legend(['Residual, rms = ',num2str(100*rms((s2-F_forward(coeff,input_nlinfit))./bsln),3),'%']);set(hleg,'box','off','location','north')
%     export_fig(['../figures/CIA_absco_',which_CIA,'_',num2str(mean([vStart,vEnd])),'.pdf'],'-q150')
    %
    
    Results{iwin,1} = coeff;
    Results{iwin,2} = [w2,s2,R]';


%%

which_CIA = 'HITRAN';
% which_CIA = 'GFIT';

P_surface = 1020;
clc;close all
input_nlinfit = [];
input_nlinfit.ils = ilsy;
optional_fields = {'shift_gas','shift_sun','shift_together','scaling','tilt','zlo'};
optional_fields_ap = [0,           0,          0,               1,       0,     0];
% the following index combined should cover
% 1:length(optional_fields)
input_fields_index = [ 3];% index as input
fit_fields_index = [1 2 4 5 6];% index included in fitting
input_fields = optional_fields(input_fields_index);
input_value = optional_fields_ap(input_fields_index);
for i = 1:length(input_value)
    input_nlinfit.(input_fields{i}) = input_value(i);
end
if ~exist('ss','var')
    ss = cell(length(window_list),1);
end

CC = lines(5);
window_list(2).wRange = [7760, 8035];
window_list(1).wRange = [9200, 9600];

addpath('~/matlab functions/export_fig/')
for iwin = 1:1%length(window_list)
    target_gas = window_list(iwin).target_gas;
    windowID = window_list(iwin).windowID;
    common_grid = window_list(iwin).common_grid;
    
    vStart = window_list(iwin).vRange(1);
    vEnd = window_list(iwin).vRange(2);
    
    wStart = window_list(iwin).wRange(1);
    wEnd = window_list(iwin).wRange(2);
    
    interval = calmat.x >= wStart & calmat.x <= wEnd;
    
    w2 = calmat.x(interval);
    s2 = double(calmat.spec(interval));
    w2 = w2(:);s2 = s2(:)/max(s2);
    
    if isempty(ss{iwin})
        int = llwaven >= vStart & llwaven <= vEnd;
        llwave = llwaven(int);lltrans = lltranswaven(int);
        ss{iwin} = F_conv_interp(llwave,lltrans,...
            2.5*common_grid_resolution,common_grid);
    end
    input_nlinfit.ss = ss{iwin};
    
    % this is very similar to mol_for_fit, except O2 CIA
    fit_spec_names = fieldnames(window_list(iwin).tau_struct);
    n_fit_spec = length(fit_spec_names);
    input_nlinfit.s1 = nan(n_fit_spec,length(common_grid));
    for i_fit_spec = 1:n_fit_spec
        tmp_tau_prof = window_list(iwin).tau_struct.(fit_spec_names{i_fit_spec}).Tau_sum;
        tmp_tau_sfc = window_list(iwin).tau_struct.(fit_spec_names{i_fit_spec}).surface_layer_sigma/1e4... % absorption cross section in m2
            *(P_surface-window_list(iwin).tau_struct.(fit_spec_names{i_fit_spec}).surface_layer_top_P)*100 ... % surface layer thickness in Pa
            /(0.029/6.02e23*9.8) ... % air molecular weight and g
            *window_list(iwin).tau_struct.(fit_spec_names{i_fit_spec}).surface_layer_VMR;
        input_nlinfit.s1(i_fit_spec,:) = double(tmp_tau_prof+tmp_tau_sfc);
    end
    
    n_fit_spec = n_fit_spec+1;
    tau_CIA = double(window_list(iwin).tau_struct.O2.(['CIA_',which_CIA]));
    input_nlinfit.s1 = cat(1,input_nlinfit.s1,tau_CIA);
    
    input_nlinfit.ils = ilsy;
    input_nlinfit.w2 = w2;
    input_nlinfit.w1 = common_grid;
    
    %             input_nlinfit.shift_gas = 0;
    %             % input_nlinfit.shift_sun = 0;
                 input_nlinfit.shift_together = 0;
    %             % input_nlinfit.scaling = 1;
    %             % input_nlinfit.tilt = 0;
    %             % input_nlinfit.zlo = 0;
    coeff0 = [optional_fields_ap(fit_fields_index) ones(1,n_fit_spec)*calmat.am];
    [coeff, R] = nlinfit(input_nlinfit,s2,@F_forward,coeff0);
    %%
%     close
    figure('unit','inch','color','w','position',[5 1 7 7])
    ax1 = subplot(2,1,1);
    bsln = polyval([coeff(4), coeff(3)],w2-mean(w2));
    h = plot(w2,s2./bsln-coeff(5),'k',w2,F_forward(coeff,input_nlinfit)./bsln-coeff(5),'linewidth',.5);
    set(h(2),'color',CC(3,:))
    hold on
%     tmp = conv(input_nlinfit.ss,ilsy/sum(ilsy),'same');
%     tmp = interp1(input_nlinfit.w1,tmp,w2,'linear','extrap');
%     plot(w2,tmp)
    for is1 = size(input_nlinfit.s1,1)
        tmp = conv(exp(-input_nlinfit.s1(is1,:)*coeff(5+is1)),ilsy/sum(ilsy),'same');
    tmp = interp1(input_nlinfit.w1,tmp,w2,'linear','extrap');
    plot(w2,tmp,'linewidth',1.5)
    end
    set(gca,'xticklabel',[])
    title(which_CIA)
    ylim([-0.05 1.05])
    hleg = legend('Observation','Fitting','CIA transmission');set(hleg,'box','off')
    if iwin == 1;set(hleg,'location','east');end
    
    ax2 = subplot(2,1,2);
    h = plot(w2,(s2-F_forward(coeff,input_nlinfit))./bsln,'k','linewidth',1.);
    pos = get(ax1,'position');
    set(ax1,'linewidth',1,'position',[pos(1) pos(2)-0.2 pos(3) pos(4)+0.2],...
        'xlim',[wStart wEnd])
    
    pos = get(ax2,'position');
    set(ax2,'linewidth',1,'position',[pos(1) pos(2) pos(3) pos(4)-0.15],...
        'xlim',[wStart wEnd])
    xlabel('Wavenumber')
    ylim([-0.05 0.05])
    hleg = legend(['Residual, rms = ',num2str(100*rms((s2-F_forward(coeff,input_nlinfit))./bsln),3),'%']);set(hleg,'box','off','location','north')
%     export_fig(['../figures/CIA_absco_',which_CIA,'_',num2str(mean([vStart,vEnd])),'.pdf'],'-q150')
    %%
    
    Results{iwin,1} = coeff;
    Results{iwin,2} = [w2,s2,R]';
end
%%
close all
subplot(2,1,1)
plot(w2,s2,'k',w2,F_forward(coeff,input_nlinfit),'linewidth',.5);
xlim([9390 9395])
subplot(2,1,2)
plot(w2,s2-F_forward(coeff,input_nlinfit),'linewidth',.5);
xlim([9390 9395])
%%
iwin = 1;
Xlim = [9350 9450];
close all
semilogy(window_list(iwin).common_grid,window_list(iwin).tau_struct.O2.Tau_sum,...
    window_list(iwin).common_grid,window_list(iwin).tau_struct.H2O.Tau_sum)
xlim(Xlim)