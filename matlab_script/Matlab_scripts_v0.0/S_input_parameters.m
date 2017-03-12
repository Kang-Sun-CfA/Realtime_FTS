%%%%%% Define input parameters
% define when the forecast was run, used for both CFS and NAM
year_run = '2016';
month_run = '08';
day_run = '04';
hour_run = '00';
% define the forecast time, should be ~noon of the measurement day
year_forecast = '2016';
month_forecast = '07';
day_forecast = '15';
hour_forecast = '18';
% define hour ahead in the forecast, an array, only used in NAM
hour_ahead = 10:24;

target_lat = 42.36;
target_lon = -71.05;

%%%%%% Define the directories
% home directory
homepath = '/home/kangsun/FTS/';
% where the matlab functions are, most importantly, the NCTOOLBOX
functionpath = '/home/kangsun/matlab functions';
% where to save the downloaded Climate Forecast data
CFSdatapath = '/data/tempo1/Shared/kangsun/FTS_data/CFS_data/';
% where to save the downloaded NAM data, prefer to using NAM
% note, the slash in the end is critical
NAMdatapath = ['/data/tempo1/Shared/kangsun/FTS_data/NAM_data/',...
year_run,month_run,day_run,'/'];
% HITRAN molecule parameter file
molparampath = '/home/kangsun/Courses/molparam.txt';
% Dir to save .mat intermediate data files
datadir = '/data/tempo1/Shared/kangsun/FTS_data/';
% path to the HITRAN spectroscopic parameters, saved using
% S_save_hitran_mat.m and F_import_par.m
hitranpath = [datadir,'FTS_5000to12000','.mat'];
% Directory to EM27 spectra
spectrapath = '/data/tempo1/Shared/kangsun/FTS_data/EM27_data/Spectra/';
% path to saved ILS
ILSfile = [datadir,'Harvard_ILS/45_20160510_1421_275.0/ergs/ilsre.dat'];

%%%%%% Define retrieval windows
field1 = 'target_gas';
field2 = 'windowID';
value1 = {'O2','CO2','CO2','CH4','CH4'};
value2 = {1,1,2,1,2};
window_list = struct(field1,value1,field2,value2);


