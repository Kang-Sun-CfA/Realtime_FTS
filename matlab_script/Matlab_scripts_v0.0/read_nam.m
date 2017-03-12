cd('D:\Research_CfA\FTS\')
addpath(genpath('C:\Users\Kang Sun\Dropbox\matlab functions'))
setup_nctoolbox
%%
% define when the forecast was run
year_run = '2016';
month_run = '07';
day_run = '26';
hour_run = '18';
% define the forecast time, should be ~noon of the measurement day
hour_ahead = [18 21 24 27];
year_forecast = '2016';
month_forecast = '07';
day_forecast = '15';
hour_forecast = '18';
Weatherdatapath = 'D:\Research_CfA\FTS\';

% url = ['www.ftp.ncep.noaa.gov/pub/data/nccf/com/nam/prod/nam.',...
%     year_run,month_run,day_run,'/',webfilename];
% outfilename = websave([Weatherdatapath,webfilename],url);
target_lat = 42.36;
target_lon = -71.05;
%%
f = ftp('ftp.ncep.noaa.gov:21');
for ihr = 1:length(hour_ahead)
webfilename = ['nam.t',hour_run,'z.awip12',num2str(hour_ahead(ihr)),...
    '.tm00.grib2'];
disp(['Downloading ',webfilename,'...'])
tic
mget(f,['pub/data/nccf/com/nam/prod/nam.',...
    year_run,month_run,day_run,'/',webfilename],...
    Weatherdatapath);
toc
end
close(f)

%'See http://www.mathworks.com/matlabcentral/fileexchange/6626-passive-mode-ftp-in-matlab'
%%
for ihr = 1:length(hour_ahead)
webfilename = ['nam.t',hour_run,'z.awip12',num2str(hour_ahead(ihr)),...
    '.tm00.grib2'];
filename = [Weatherdatapath,'pub\data\nccf\com\nam\prod\nam.',...
    year_run,month_run,day_run,'\',webfilename];
nco = ncgeodataset(filename);

end
%%
filename = 'nam.t00z.conusnest.hiresf00.tm00.grib2';
nco = ncgeodataset(filename);
nco.variables
namData = F_read_nam(filename,target_lat,target_lon);
%ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/nam/prod/nam.20160726/
%%
param = 'x';
x = nco{param}(:);
param = 'y';
y = nco{param}(:);
param = 'Pressure_surface';
Pressure_surface = nco{param}(:);
temp = nco{'LambertConformal_Projection'};temp.attributes

param = 'Temperature_surface';
Temperature_surface = nco{param}(:);
%%
param = 'Temperature_isobaric';
temp = nco.geovariable(param);
Temperature_isobaric = temp(1,:,5:9,:);
%%
close all
ax1 = subplot(1,2,1);
h = pcolor(namData.x(int_x),namData.y(int_y),double(squeeze(namData.Temperature_isobaric(1,:,:))));
set(h,'edgecolor','none')
ax2 = subplot(1,2,2);
h = pcolor(namData.x,namData.y,double((squeeze(namData.Temperature_surface(1,:,:)))));
set(h,'edgecolor','none')
set(ax2,'clim',get(ax1,'clim'),'xlim',get(ax1,'xlim'),'ylim',...
    get(ax1,'ylim'))

%%
clc
mstruct = defaultm('lambertstd');
mstruct.origin = [25 265 0];
mstruct.mapparallels = 25;
mstruct.nparallels = 1;
% scale factor = earth radius, in km
mstruct.scalefactor = 6371229/1e3;
mstruct.falseeasting = 0;
mstruct.falsenorthing = 0;
mstruct = defaultm(mstruct);
% coordinate of Boston
target_lon = -71.06;
target_lat = 42.36;
[target_x,target_y] = mfwdtran(mstruct,target_lat,target_lon);
NSA = shaperead('landareas.shp', 'UseGeoCoords', true,...
    'Selector',{@(name)strcmpi(name,'North and South America'),'Name'});
[NSA_x,NSA_y] = mfwdtran(mstruct,NSA.Lat,NSA.Lon);
%%
close all;figure('color','w','unit','inch','position',[0 0.5 7.5 6])
opengl software
h = pcolor(namData.x,namData.y,double(squeeze(namData.Temperature_surface)));
set(h,'edgecolor','none')
hold on;plot(target_x,target_y,'-pr','markersize',12,'linewidth',1);
plot(NSA_x,NSA_y,'w','linewidth',2)
hold off;text(target_x-900,target_y+200,'Boston, where I want T/P')
xlabel('Distance [km]');ylabel('Distance [km]')
hc = colorbar;set(get(hc,'ylabel'),'string','Surface T [K]')
title(['NAM forecast made at ',month_run,day_run,hour_run])