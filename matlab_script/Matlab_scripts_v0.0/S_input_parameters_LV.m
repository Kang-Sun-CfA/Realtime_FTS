%%%%%% Define input parameters
% define when the forecast was run
year_run = '2016';
month_run = '08';
day_run = '04';
hour_run = '00';

% define hour ahead in the forecast, an array
hour_ahead = 10:24;

target_lat = 42.36;
target_lon = -71.05;

%%%%%% Define the directories
% main data working directory
datadir = 'D:\Research_CfA\FTS\';

% where to save the downloaded NAM data, prefer to using NAM
% note, the slash in the end is critical
NAMdatapath = [datadir,'NAM_data\',year_run,month_run,day_run,'\'];

% path to the HITRAN spectroscopic parameters, saved using
% S_save_hitran_mat.m and F_import_par.m
% hitranpath = [datadir,'FTS_5000to12000','.mat'];
% Directory to EM27 spectra
% spectrapath = '/data/tempo1/Shared/kangsun/FTS_data/EM27_data/Spectra/';
% path to saved ILS
% ILSfile = [datadir,'Harvard_ILS/45_20160510_1421_275.0/ergs/ilsre.dat'];

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
