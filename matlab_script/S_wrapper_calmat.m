%% IMPORTANT INPUTS
%%%%%% Define retrieval windows
field1 = 'target_gas';
field2 = 'windowID';
value1 = {'O2','CO2','CO2','CH4','CH4'};
value2 = {1,1,2,1,2};
window_list = struct(field1,value1,field2,value2);

%%%%%% Define max wait time for a new file, in seconds
maxquerytime = 5;

%%%%%% Define the resolution of forward model, should be consistant with HITRAN
common_grid_resolution = 0.05;

%%%%%% Use input ILS calculated in Labview (1) or every file (0)
use_input_ILS = 1;

%%%%%% Half length of the ILS, the longer the more accurate and slower
ilsextent = 10; % wavenumber

%%%%%% system separator, "/" for mac and linux, "\" for windows
sfs = filesep;

%%%%%% rate of downsample. you don't need smoothed ifg at original rate
ifgs_step = 500;
%% INPUTS FOR CALMAT
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
%% INPUT FROM LABVIEW WRAPPER USED FOR TESTING, SHOULD BE COMMENTED OUT
date_obs = '2016,06,24';
observation = 'hb';
fileN = 1138;
last_fileN = -1;
spec_path = 'C:\Users\ksun\Documents\GitHub\Realtime_FTS\spectra\160624';
do_fit = false;
%% FIND THE CORRECT FILENAME
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
ZPDRight = calmat.ZPDRight;
ZPDLeft = calmat.ZPDLeft;
ifg_smooth = calmat.ifg_smooth(1:ifgs_step:end);
%% SEPARATE INTO FITTING WINDOWS AND OPTIONALLY DO FITTING
Results = cell(length(window_list),2);
if ~exist('do_fit','var')
do_fit = true;
end
for iwin = 1:length(window_list)
    target_gas = window_list(iwin).target_gas;
    windowID = window_list(iwin).windowID;
    [vStart, vEnd, wStart, wEnd, mol_for_fit,mol_for_spec] = ...
        F_define_windows(target_gas,windowID);
    interval = calmat.x >= wStart & calmat.x <= wEnd;
    
    w2 = calmat.x(interval);
    s2 = double(calmat.spec(interval));
    w2 = w2(:);s2 = s2(:);
    if ~do_fit;R = nan*s2;coeff = nan;end;
    Results{iwin,1} = coeff;
    Results{iwin,2} = [w2,s2,R]';
end