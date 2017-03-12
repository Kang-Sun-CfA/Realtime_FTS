function namData = F_read_nam(filename,target_lat,target_lon)
% this function reads NAM data given the filename, outputs temperautre,
% pressure, height profiles at the target coordinate
% written by Kang Sun on 2016/7/29

nco = ncgeodataset(filename);

param = 'x';
namData.x = nco{param}(:);

param = 'y';
namData.y = nco{param}(:);

mstruct = defaultm('lambertstd');
mstruct.origin = [25 265 0];
mstruct.mapparallels = 25;
mstruct.nparallels = 1;
mstruct.scalefactor = 6371229/1e3;
mstruct.falseeasting = 0;
mstruct.falsenorthing = 0;
mstruct = defaultm(mstruct);

[target_x, target_y] = mfwdtran(mstruct,target_lat,target_lon);
int_x = find(namData.x > target_x-1000 & namData.x < target_x+1000);
int_y = find(namData.y > target_y-2000 & namData.y < target_y+2000);
param = 'isobaric3';
namData.isobaric = nco{param}(:);
param = 'Temperature_isobaric';
temp = nco.geovariable(param);
namData.Temperature_isobaric = squeeze(temp(1,:,int_y,int_x));
param = 'Geopotential_height_isobaric';
temp = nco.geovariable(param);
namData.Geopotential_height_isobaric = squeeze(temp(1,:,int_y,int_x));
param = 'Specific_humidity_isobaric';
temp = nco.geovariable(param);
namData.Specific_humidity_isobaric = squeeze(temp(1,:,int_y,int_x));

namData.namT = interptoxy(namData.Temperature_isobaric,namData.x(int_x),namData.y(int_y),...
target_x,target_y,'linear');
namData.namH = interptoxy(namData.Geopotential_height_isobaric,namData.x(int_x),namData.y(int_y),...
target_x,target_y,'linear');
namData.namH2O = interptoxy(namData.Specific_humidity_isobaric,namData.x(int_x),namData.y(int_y),...
target_x,target_y,'linear')/0.622;

param = 'Pressure_surface';
namData.Pressure_surface = nco{param}(:);
param = 'Temperature_surface';
namData.Temperature_surface = nco{param}(:);
%param = 'Specific_humidity_isobaric';
%namData.Specific_humidity_isobaric = nco{param}(:);
