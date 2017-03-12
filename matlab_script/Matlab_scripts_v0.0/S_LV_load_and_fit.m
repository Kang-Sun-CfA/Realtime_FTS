% functionpath = 'D:\Research_CfA\FTS\Matlab_scripts';
% spec_path = 'D:\Research_CfA\FTS\ha_20150807';
% ifg_path = 'D:\Research_CfA\FTS\ha_20150807_ifg';
% observation = 'ha';
% fileN = 100;
clc
S_input_parameters_LV;
temp = textscan(date_obs,'%s%s%s','delimiter',',');

obs_year = temp{1}{1};
obs_month = temp{2}{1};
obs_day = temp{3}{1};
% it is tricky how the obus files are named, and whether ifg are saved together with spectra
if str2double(obs_year) >= 2016
SBstring = 's0e00a.';
else
SBstring = '.s0e00a.';
end
datenum([obs_year,'-',obs_month,'-',obs_day])

filename = [spec_path,'\',observation,[obs_year,obs_month,obs_day],SBstring,sprintf('%0*d',4,fileN)];

if last_fileN >= 0 && fileN > last_fileN
nextdaystr = datestr(datenum([obs_year,'-',obs_month,'-',obs_day])+1,'yyyymmdd');

filename = [spec_path,'\',observation,nextdaystr,SBstring,sprintf('%0*d',4,fileN-last_fileN)];

ifgname = [ifg_path,'\',observation,nextdaystr,'.ifg.',sprintf('%0*d',4,fileN-last_fileN)];
end

if str2double(obs_year) >= 2016
    ifgname = filename;
else
ifgname = [ifg_path,'\',observation,[obs_year,obs_month,obs_day],'.ifg.',sprintf('%0*d',4,fileN)];
end

tic
while ~exist(filename,'file')  
    pause(1)
    querytime = toc;
    if querytime > maxquerytime
        error(['Cannot find OPUS file #',num2str(fileN),' in ',num2str(maxquerytime),' s!'])
    end
end
[spectrum,v, params] = ImportOpus(filename,'SampleSpectrum');
TSC = params.Instrument.TSC;
DAT = params.SampleDataSpectrum.DAT(1:end-2);
TIM = params.SampleDataSpectrum.TIM(1:end-4);

if ~exist(ifgname,'file')
ifg = 0;t = 0;
else
[ifg,t] = ImportOpus(ifgname,'SampleInterferogram');
ifg = double(ifg(:)');
end
spectrum = double(spectrum(:)');
obs_x = v;
obs_y = spectrum;
fclose('all');
Results = cell(length(window_list),2);
if ~exist('do_fit','var')
do_fit = true;
end
if prctile(spectrum(v > 5900 & v < 6500),90) < 0.06
    do_fit = false;
end
load([NAMdatapath,'fcst_prof',year_run,month_run,day_run,hour_run,'.mat'])
load([datadir, 'FTS_solar.mat'],'llwaven','lltranswaven')
ihr = 1;
Z_level = namH_mat(:,ihr);
dZ = diff(Z_level);
P_level = namP_mat(:,ihr)/100; % pressure in hPa
P_layer = P_level(1:end-1)+diff(P_level)/2;
nlayer = length(P_layer);
T_level = namT_mat(:,ihr);
T_layer = interp1(P_level,T_level,P_layer);
Z_layer = interp1(P_level,Z_level,P_layer);
savepath = [datadir,'Tau_data\',year_run,...
    month_run,day_run,'\UTC',hour_run,'f',num2str(hour_ahead(ihr)),'\'];

for iwin = 1:length(window_list)
    target_gas = window_list(iwin).target_gas;
    windowID = window_list(iwin).windowID;
    [vStart, vEnd, wStart, wEnd, mol_for_fit,mol_for_spec] = ...
        F_define_windows(target_gas,windowID);
    tempS = load([savepath,'common_grid_',num2str(vStart),...
        '_',num2str(vEnd),'.mat'],'common_grid');
    common_grid = tempS.common_grid;
    if abs(median(diff(common_grid))-common_grid_resolution) > 1e-10
        error('Common_grid_resolution is not consistent! Check S_input_parameters_LV.m')
    end
    % common_grid_resolution = median(diff(common_grid));
    dtau_mol = nan(length(mol_for_fit),length(common_grid),nlayer);
    % clear the fitting input from the previous window
    input = struct;
    for imol = 1:length(mol_for_fit)
        taumol = mol_for_fit{imol};
        if strcmpi(taumol,'HDO')
            tempS = load([savepath,'tau_data_','H2O','_',num2str(vStart),...
                '_',num2str(vEnd),'.mat'],'output_LBL');
            output_LBL = tempS.output_LBL;
            dtau_mol(imol,:,:) = output_LBL.dtau_layer(2,:,:);
        elseif strcmpi(taumol,'O2')
            tempS = load([savepath,'tau_data_',taumol,'_',num2str(vStart),...
                '_',num2str(vEnd),'.mat'],'output_LBL');
            output_LBL = tempS.output_LBL;
            dtau_mol(imol,:,:) = output_LBL.dtau_layer(1,:,:);
            % static input to forward function
            input.s_CIA = double(sum(output_LBL.tau_cntm,1));
        else
            tempS = load([savepath,'tau_data_',taumol,'_',num2str(vStart),...
                '_',num2str(vEnd),'.mat'],'output_LBL');
            output_LBL = tempS.output_LBL;
            dtau_mol(imol,:,:) = output_LBL.dtau_layer(1,:,:);
        end
    end
    int = llwaven >= vStart & llwaven <= vEnd;
    llwave = llwaven(int);lltrans = lltranswaven(int);
    int = obs_x >= wStart & obs_x <= wEnd;
    w2 = double(obs_x(int))';
    s2 = double(obs_y(int))';
    
    % use ideal, boxcar ILS.
    if use_input_ILS
        ilsy = ILS_input(iwin,:);
    else
        vCenter = mean([vStart vEnd]);
        %     ilsextent = 10; % wavenumber
        ilsx = -ilsextent:common_grid_resolution:ilsextent;
        OPD = 1.8; % maximum opd in cm
        FOV = 2.96e-3; % field of view in rad
        ils_sync = 2*OPD*sinc(2*ilsx*OPD);
        if sum(ilsx >= -.5*vCenter*FOV^2 & ilsx <= 0) > 1
            ils_rect = zeros(size(ilsx));
            ils_rect(ilsx >= -.5*vCenter*FOV^2 & ilsx <= 0) = 2/(vCenter*FOV^2);
            ilsy = conv(ils_sync,ils_rect,'same')*common_grid_resolution;
        else
            ilsy = ils_sync;
        end
        modulation_loss = 0.025;
        ils_align = sinc(ilsx*OPD/modulation_loss).^2*OPD/abs(modulation_loss);
        
        ilsy = conv(ilsy,ils_align,'same')*common_grid_resolution;
    end
    % static input to forward function
    input.ils = ilsy;
    input.w2 = w2;
    input.w1 = common_grid;
    input.ss = ...
        F_conv_interp(llwave,lltrans,...
        2.5*common_grid_resolution,common_grid);
    input.s1 = squeeze(sum(dtau_mol,3));
    if strcmpi(target_gas,'O2')
        % tilt is important, ZLO is not
        input.varlist = {'cntm_level','cntm_tilt',...
            'gas_shift',...
            'solar_scaling','CIA_scaling',...
            'gas_scaling'};
        coeff0 = [prctile(s2,90) 0,...
            0,...
            1, 0.17,...
            ones(1,length(mol_for_fit))];
    else
        input.varlist = {'cntm_level','cntm_tilt',...
            'gas_shift',...
            'solar_scaling',...
            'gas_scaling'};
        coeff0 = [prctile(s2,90) 0,...
            0,...
            1,...
            ones(1,length(mol_for_fit))];
    end
    % all possible state vector:
    % 1 continuum level, 2 continuum tilt
    % 3 gas wave number shift, 4 solar wavenumber shift, 5 ZLO
    % 6 solar scaling, 7 CIA scaling
    % 8+ add target gas scaling
    % coeff0 = [prctile(s2,90) 0,...
    %  0 0 1,...
    %  1, 1,...
    %  ones(1,length(mol_for_fit)) 0.17];
    if do_fit
try
        [coeff, R] = nlinfit(input,s2,@F_forward_FTS,coeff0);
catch
coeff = coeff0*0;
R = nan*s2;
end
    else
        coeff = coeff0*0;
        R = nan*s2;
    end
    Results{iwin,1} = coeff;
    % R is defined as: R = s2-F_forward_FTS(coeff,input)
    Results{iwin,2} = [w2,s2,R]';
end
