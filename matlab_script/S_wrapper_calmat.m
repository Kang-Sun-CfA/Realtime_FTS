% system separator, "/" for mac and linux, "\" for windows
sfs = filesep;
%% INPUT FROM LABVIEW WRAPPER USED FOR TESTING, SHOULD BE COMMENTED OUT
% date_obs = '2016,06,24';
% observation = 'hb';
% fileN = 1856;
% last_fileN = -1;
% spec_path = ['..',sfs,'spectra',sfs,'160624'];
% do_fit = true;

ifamf = true;
which_CIA = 'GFIT';% GFIT, SAO, or Mate
%% IMPORTANT INPUTS
% currently do not support updating optical depth on the fly, use presaved
update_profile = false;
if ~exist('do_fit','var')
    do_fit = true;
end
% Load fitting window information and pre-calculated optical depth
profile_date = '0000';% mmdd, 0000 for US standard atm
profile_hour = '0000';% hhmm, 0000 for US standard atm
window_list_fn = ['..',sfs,'spectroscopy',sfs,'window_list_',...
    profile_date,'_',profile_hour,'.mat'];

if do_fit && update_profile
    inp_tau = [];
    
    %!!!!!!!!!!!!!!!! local dir saving ABSCO tables !!!!!!!!!!!!!!!!!!!!!!
    inp_tau.absco_dir = '/data/tempo1/Shared/kangsun/FTS_data/HITRAN/';
    %!!!!!!!!!!!!!!!! YOU HAVE TO SPECIFY THIS !!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    inp_tau.CIA_dir = ['..',sfs,'spectroscopy',sfs];
    inp_tau.which_CIA = {'sao','gfit','mate'};
    inp_tau.profile_fn = ['..',sfs,'profiles',sfs,'profiles_',profile_date,...
        '_',profile_hour,'.dat'];
    inp_tau.lowest_possible_Psurf = 980;
    inp_tau.nsublayer = 3;
    inp_tau.common_grid_resolution = 0.05;
    inp_tau.surface_layer_P = 990;
    
    window_list = F_absco2tau(inp_tau);
    
    save(window_list_fn,'window_list')
end
% Define retrieval windows
if exist(window_list_fn,'file')
    load(window_list_fn)
    common_grid_resolution = mean(diff(window_list(1).common_grid));
else
    warning('Window list file does not exist!')
    field1 = 'target_gas';
    field2 = 'windowID';
    value1 = {'O2','CO2','CO2','CH4','CH4'};
    value2 = {1,1,2,1,2};
    window_list = struct(field1,value1,field2,value2);
    common_grid_resolution = 0.05;
end

% Define max wait time for a new file, in seconds
maxquerytime = 5;

% which apodization.
%        < 0: use the OPUS spectrum, no fft on the ifg
%          1: boxcar
%          2: triag
%          3: Happ-Genzel
%        4/5: Blackmann-Harris 3-term/4-term
%      6/7/8: Norton-Beer weak/medium/strong
% apokind_all = 7;

% rate of downsample. you don't need smoothed ifg at original rate
ifgs_step = 500;
%% WORK OUT AN ILS, if necessary
if ~exist('ilsy','var')
    % Half length of the ILS, the longer the more accurate and slower
    ilsextent = 10; % wavenumber
    % find an OPUS file with the same ILS as the measurements
    ils_fn = ['..',sfs,'spectra',sfs,'160624',sfs,'hb20160624s0e00a.1138'];
    input_ils = [];
    input_ils.filename = ils_fn;
    input_ils.sim_wn = 1.579792480468750e+04/2;
    input_ils.phase_error = 0*5.412e-1;
    input_ils.apokind = apokind_all; % should be consistent with apokind in config
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
end
%% INPUTS FOR CALMAT
% initialize the input
config = [];

config.apokind = apokind_all;

% lower and upper limit of output spectrum, in cm-1
config.lowerlim = 5000;
config.upperlim = 10000;

% display processing message or not. I'd like it to be off
config.display = false;

config.ifamf = ifamf;
if do_fit
    config.ifamf = true;
end

config.loclat = 39.767;
config.loclon = -86.15;
config.localt = 0.22;

% time difference compared to UTC in hour
config.timezone = -5;
%% FIND THE CORRECT FILENAME
if last_fileN < 0 || (last_fileN >= 0 && fileN <= last_fileN)
    temp = textscan(date_obs,'%s%s%s','delimiter',',');
    
    obs_year = temp{1}{1};
    obs_month = temp{2}{1};
    obs_day = temp{3}{1};
    % it is tricky how the obus files are named, and whether ifg are saved together with spectra
    % this is an ad hoc solution
    if str2double(obs_year) >= 2016
        SBstring = 's0e00a.';
    else
        SBstring = '.s0e00a.';
    end
    % datenum([obs_year,'-',obs_month,'-',obs_day])
    filename = [spec_path,sfs,observation,[obs_year,obs_month,obs_day],...
        SBstring,sprintf('%0*d',4,fileN)];
    
    % if the filename does not exist, wait for maxquerytime (s)
    tic
    while ~exist(filename,'file')
        pause(1)
        querytime = toc;
        if querytime > maxquerytime
            error(['Cannot find OPUS file #',num2str(fileN),' in ',num2str(maxquerytime),' s!'])
        end
    end
    %% F. F. T!
    calmat = F_calmat(filename,config);
    %% PICK OUT SEVERAL KEY OUTPUT VARIABLES
    ifgQuality = calmat.ifgQuality;
    percentDeviation = calmat.percentDeviation;
    ZPDRight = calmat.ZPDRight';
    ZPDLeft = calmat.ZPDLeft';
    ifg_smooth = calmat.ifg_smooth(1:ifgs_step:end)';
    AMF = calmat.am;
    nburstfwd = calmat.nburstfwd;
    nburstbwd = calmat.nburstbwd;
    scanT = calmat.scanT;
    %% SEPARATE INTO FITTING WINDOWS AND OPTIONALLY DO FITTING
    Results = cell(length(window_list),2);
    % prepare data base if do fitting
    if do_fit
        load(['..',sfs,'spectroscopy',sfs,'FTS_solar.mat'],'llwaven','lltranswaven')
        input_nlinfit = [];
        input_nlinfit.ils = ilsy;
        optional_fields = {'shift_gas','shift_sun','shift_together','scaling','tilt','zlo'};
        optional_fields_ap = [0,           0,          0,               1,       0,     0];
        % the following index combined should cover
        % 1:length(optional_fields)
        input_fields_index = [1 3];% index as input
        fit_fields_index = [2 4 5 6];% index included in fitting
        input_fields = optional_fields(input_fields_index);
        input_value = optional_fields_ap(input_fields_index);
        for i = 1:length(input_value)
            input_nlinfit.(input_fields{i}) = input_value(i);
        end
        if ~exist('ss','var')
            ss = cell(length(window_list),1);
        end
    end
    for iwin = 1:length(window_list)
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
        if ~do_fit
            R = nan*s2;coeff = nan;
        else
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
                input_nlinfit.s1(i_fit_spec,:) = ...
                    double(window_list(iwin).tau_struct.(fit_spec_names{i_fit_spec}).Tau_sum);
            end
            % add O2 CIA for the O2 fitting window
            if iwin == 1
                n_fit_spec = n_fit_spec+1;
                tau_CIA = double(window_list(iwin).tau_struct.O2.(['CIA_',which_CIA]));
                input_nlinfit.s1 = cat(1,input_nlinfit.s1,tau_CIA);
            end
            input_nlinfit.ils = ilsy;
            input_nlinfit.w2 = w2;
            input_nlinfit.w1 = common_grid;
            
%             input_nlinfit.shift_gas = 0;
%             % input_nlinfit.shift_sun = 0;
%             input_nlinfit.shift_together = 0;
%             % input_nlinfit.scaling = 1;
%             % input_nlinfit.tilt = 0;
%             % input_nlinfit.zlo = 0;
            coeff0 = [optional_fields_ap(fit_fields_index) ones(1,n_fit_spec)*AMF];
            [coeff, R] = nlinfit(input_nlinfit,s2,@F_forward,coeff0);
%             figure
%             plot(w2,s2,'k',w2,F_forward(coeff0,input_nlinfit),w2,F_forward(coeff,input_nlinfit))

        end
        Results{iwin,1} = coeff;
        Results{iwin,2} = [w2,s2,R]';
    end
end