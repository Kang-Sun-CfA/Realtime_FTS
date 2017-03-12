clc

% load input parameters
S_input_parameters;

cd(homepath)
addpath(genpath(functionpath))

webfilename = ['pgbf',year_forecast,month_forecast,day_forecast,hour_forecast,...
    '.01.',year_run,month_run,day_run,hour_run,'.grb2'];
url = ['http://nomads.ncdc.noaa.gov/thredds/fileServer/modeldata/cfsv2_forecast_6-hourly_9mon_pgbf/'...
    year_run,'/',year_run,month_run,'/',year_run,month_run,day_run,'/',...
    year_run,month_run,day_run,hour_run,...
    '/',webfilename];

outfilename = websave([CFSdatapath,webfilename],url);

% Get the forecast temperature and pressure profiles!
try
    [cfsP,cfsH,cfsT,cfsH2O] = F_read_cfs(outfilename,target_lat,target_lon);
catch
    setup_nctoolbox
    [cfsP,cfsH,cfsT,cfsH2O] = F_read_cfs(outfilename,target_lat,target_lon);
end

% assume geometric height is geopotential height, Is that OK?
Z_level = cfsH; 
dZ = diff(Z_level);
P_level = cfsP/100; % pressure in hPa
P_layer = P_level(1:end-1)+diff(P_level)/2;
T_level = cfsT;
T_layer = interp1(P_level,T_level,P_layer);
Z_layer = interp1(P_level,Z_level,P_layer);
%% Hack GFIT vmr profiles
fid = fopen('/data/tempo1/Shared/kangsun/ggg_location/vmrs/gnd/gnd_summer.vmr');
C = textscan(fid,['%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f',...
'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f'],...
'headerlines',5);
C = cell2mat(C);
fclose(fid);
profileinput.C = C;
profileinput.Z_layer = Z_layer;
profileinput.cfsH2O = cfsH2O;
%%
target_gas = 'O2';
windowID = 1;
%%
target_gas = 'H2O';
windowID = 1;
%%
target_gas = 'CO2';
windowID = 1;
%%
clc
[vStart, vEnd, wStart, wEnd, molnameinput,mol_for_spec] = ...
    F_define_windows(target_gas,windowID);
for imol = 1:length(mol_for_spec)
    taumol = mol_for_spec{imol};
    [w_data, tau_data] = ...
    F_line_by_line_FTS(taumol,F_ap_profile(taumol,profileinput),P_level,...
    P_layer,T_level,T_layer,Z_level,Z_layer,vStart,vEnd,molparampath,datadir);
end
% %% the O2 fitting window
% vStart = 7700; vEnd = 8050;
% % O2, lines
% taumol = 'O2';
% [w_data, tau_data] = ...
%     F_line_by_line_FTS(taumol,F_ap_profile(taumol,profileinput),P_level,...
%     P_layer,T_level,T_layer,Z_level,Z_layer,vStart,vEnd,molparampath,datadir);
% % O2, CIA. The O2-N2 CIA looks quite weird, not use. Maybe HITRAN's problem
% taumol = 'O2_cntm';
% [w_data, tau_data] = ...
%     F_line_by_line_FTS(taumol,F_ap_profile('O2',profileinput),P_level,...
%     P_layer,T_level,T_layer,Z_level,Z_layer,vStart,vEnd,molparampath,datadir);
% % H2O
% taumol = 'H2O';
% [w_data, tau_data] = ...
%     F_line_by_line_FTS(taumol,F_ap_profile(taumol,profileinput),P_level,...
%     P_layer,T_level,T_layer,Z_level,Z_layer,vStart,vEnd,molparampath,datadir);
% % CO2
% taumol = 'CO2';
% [w_data, tau_data] = ...
%     F_line_by_line_FTS(taumol,F_ap_profile(taumol,profileinput),P_level,...
%     P_layer,T_level,T_layer,Z_level,Z_layer,vStart,vEnd,molparampath,datadir);
% %% the CO2 fitting window, # 1: CO2, H2O, HDO, CH4
% % CO2
% vStart = 6160; vEnd = 6280;
% taumol = 'CO2';
% [w_data, tau_data] = ...
%     F_line_by_line_FTS(taumol,F_ap_profile(taumol,profileinput),P_level,...
%     P_layer,T_level,T_layer,Z_level,Z_layer,vStart,vEnd,molparampath,datadir);
% % H2O, HDO is included
% taumol = 'H2O';
% [w_data, tau_data] = ...
%     F_line_by_line_FTS(taumol,F_ap_profile(taumol,profileinput),P_level,...
%     P_layer,T_level,T_layer,Z_level,Z_layer,vStart,vEnd,molparampath,datadir);
% % CH4
% taumol = 'CH4';
% [w_data, tau_data] = ...
%     F_line_by_line_FTS(taumol,F_ap_profile(taumol,profileinput),P_level,...
%     P_layer,T_level,T_layer,Z_level,Z_layer,vStart,vEnd,molparampath,datadir);
% %% the CH4 fitting window, # 1:  CH4, CO2, H2O, N2O,
% % CO2
% wStart = 5880; wEnd = 5996;
% vStart = wStart-20; vEnd = wEnd+20;
% taumol = 'CO2';
% [w_data, tau_data] = ...
%     F_line_by_line_FTS(taumol,F_ap_profile(taumol,profileinput),P_level,...
%     P_layer,T_level,T_layer,Z_level,Z_layer,vStart,vEnd,molparampath,datadir);
% taumol = 'H2O';
% [w_data, tau_data] = ...
%     F_line_by_line_FTS(taumol,F_ap_profile(taumol,profileinput),P_level,...
%     P_layer,T_level,T_layer,Z_level,Z_layer,vStart,vEnd,molparampath,datadir);
% taumol = 'CH4';
% [w_data, tau_data] = ...
%     F_line_by_line_FTS(taumol,F_ap_profile(taumol,profileinput),P_level,...
%     P_layer,T_level,T_layer,Z_level,Z_layer,vStart,vEnd,molparampath,datadir);
% taumol = 'N2O';
% [w_data, tau_data] = ...
%     F_line_by_line_FTS(taumol,F_ap_profile(taumol,profileinput),P_level,...
%     P_layer,T_level,T_layer,Z_level,Z_layer,vStart,vEnd,molparampath,datadir);
% %% the CH4 fitting window, # 2:  CH4, CO2, H2O, HDO,
% % CO2
% wStart = 6007; wEnd = 6145;
% vStart = wStart-20; vEnd = wEnd+20;
% taumol = 'CO2';
% [w_data, tau_data] = ...
%     F_line_by_line_FTS(taumol,F_ap_profile(taumol,profileinput),P_level,...
%     P_layer,T_level,T_layer,Z_level,Z_layer,vStart,vEnd,molparampath,datadir);
% taumol = 'H2O';
% [w_data, tau_data] = ...
%     F_line_by_line_FTS(taumol,F_ap_profile(taumol,profileinput),P_level,...
%     P_layer,T_level,T_layer,Z_level,Z_layer,vStart,vEnd,molparampath,datadir);
% taumol = 'CH4';
% [w_data, tau_data] = ...
%     F_line_by_line_FTS(taumol,F_ap_profile(taumol,profileinput),P_level,...
%     P_layer,T_level,T_layer,Z_level,Z_layer,vStart,vEnd,molparampath,datadir);
% %%
% close all;hold on
% for ilayer = 1:10:length(P_layer);
%     plot(w_data{ilayer},tau_data{ilayer})
%     set(gca,'yscale','log')
% end
% %%
% close all;hold on
% for ilayer = 1:10:length(P_layer);
%     plot(w_data,tau_data{ilayer})
%     set(gca,'yscale','log')
% end
% %%
% semilogx(C(:,12),C(:,1))