function [cfsP,cfsH,cfsT,cfsH2O] = F_read_cfs(filename,target_lat,target_lon)
% this function reads CFS data given the filename, outputs temperature,
% pressure, geopotential height profile at the target coordinate
% written by Kang Sun on 2016/07/15

% url = ['http://nomads.ncdc.noaa.gov/modeldata/cfsv2_forecast_mm_9mon/',...
% '2016/201607/20160714/2016071418/',...
%     'pgbf.01.2016071418.201612.avrg.grib.grb2'];
% 
% if ~exist('ncgeodataset')
%     error('You need to install NCTOOLBOX to read GRIB2 data')
% end
% clc
% % filename = ['/data/tempo1/Shared/kangsun/FTS_data/',...
% %     'pgbf.01.2016071418.201612.avrg.grib.grb2'];
% filename = ['/data/tempo1/Shared/kangsun/FTS_data/',...
%     'pgbf2016071600.01.2016071418.grb2'];

nco=ncgeodataset(filename);

% target_lat = 42.36;
% target_lon = -71.05;
if target_lon < 0
    target_lon = 360+target_lon;
end
lat=nco{'lat'}(:);
lon=nco{'lon'}(:);

param = 'Temperature_isobaric';
Temperature_isobaric = squeeze(nco{param}(:));
Temperature_isobaric = permute(Temperature_isobaric,[2,3,1]);

param = 'Specific_humidity_isobaric';
Specific_humidity_isobaric = squeeze(nco{param}(:));
Specific_humidity_isobaric = permute(Specific_humidity_isobaric,[2,3,1]);

param = 'isobaric';
isobaric = nco{param}(:);

Geopotential_height_isobaric = squeeze(nco{'Geopotential_height_isobaric'}(:));
Geopotential_height_isobaric = permute(Geopotential_height_isobaric,[2,3,1]);

[longrid, latgrid] = meshgrid(lon,lat);
nlevel = length(isobaric);
cfsP = isobaric;
cfsT = isobaric;
cfsH = isobaric;
cfsH2O = isobaric;
for ilevel = 1:nlevel
    cfsT(ilevel) = ...
        interp2(longrid,latgrid,squeeze(Temperature_isobaric(:,:,ilevel)),target_lon,target_lat);
    cfsH(ilevel) = ...
        interp2(longrid,latgrid,squeeze(Geopotential_height_isobaric(:,:,ilevel)),target_lon,target_lat);
    cfsH2O(ilevel) = ...
        interp2(longrid,latgrid,squeeze(Specific_humidity_isobaric(:,:,ilevel)),target_lon,target_lat);
end
cfsH2O = cfsH2O/0.622; % change from mass to volume ratio
% plot(cfsT,cfsH/1000,'-o')
% %%
% h = pcolor(longrid,latgrid,squeeze(Temperature_isobaric(:,:,1)));
% set(h,'edgecolor','none');colorbar
% %%
% nco.variables
% 
% ans = 
% 
%     'pressure_difference_layer_bounds'
%     'pressure_difference_layer1_bounds'
%     'sigma_layer_bounds'
%     'height_above_ground_layer_bounds'
%     'height_above_ground_layer1_bounds'
%     'time1_bounds'
%     'Temperature_maximum_wind'
%     'Temperature_tropopause'
%     'Temperature_isobaric'
%     'Temperature_altitude_above_msl'
%     'Temperature_sigma'
%     'Temperature_pressure_difference_layer'
%     'Temperature_potential_vorticity_surface'
%     'Potential_temperature_sigma'
%     'Dew-point_temperature_height_above_ground'
%     'Dew-point_temperature_pressure_difference_layer'
%     'Specific_humidity_isobaric'
%     'Specific_humidity_pressure_difference_layer'
%     'Relative_humidity_zeroDegC_isotherm'
%     'Relative_humidity_isobaric'
%     'Relative_humidity_height_above_ground'
%     'Relative_humidity_sigma'
%     'Relative_humidity_sigma_layer'
%     'Relative_humidity_pressure_difference_layer'
%     'Relative_humidity_entire_atmosphere'
%     'Relative_humidity_highest_tropospheric_freezing'
%     'Precipitable_water_pressure_difference_layer'
%     'Total_precipitation_surface_6_Hour_Accumulation'
%     'Large-scale_precipitation_non-convective_surface_6_Hour_Accumulation'
%     'Convective_precipitation_surface_6_Hour_Accumulation'
%     'Cloud_mixing_ratio_isobaric'
%     'Categorical_Rain_surface_6_Hour_Average'
%     'Categorical_Freezing_Rain_surface_6_Hour_Average'
%     'Categorical_Ice_Pellets_surface_6_Hour_Average'
%     'Categorical_Snow_surface_6_Hour_Average'
%     'u-component_of_wind_maximum_wind'
%     'u-component_of_wind_tropopause'
%     'u-component_of_wind_isobaric'
%     'u-component_of_wind_altitude_above_msl'
%     'u-component_of_wind_sigma'
%     'u-component_of_wind_pressure_difference_layer'
%     'u-component_of_wind_potential_vorticity_surface'
%     'v-component_of_wind_maximum_wind'
%     'v-component_of_wind_tropopause'
%     'v-component_of_wind_isobaric'
%     'v-component_of_wind_altitude_above_msl'
%     'v-component_of_wind_sigma'
%     'v-component_of_wind_pressure_difference_layer'
%     'v-component_of_wind_potential_vorticity_surface'
%     'Stream_function_isobaric'
%     'Velocity_potential_isobaric'
%     'Vertical_velocity_pressure_isobaric'
%     'Vertical_velocity_pressure_sigma'
%     'Absolute_vorticity_isobaric'
%     'Vertical_Speed_Shear_tropopause'
%     'Vertical_Speed_Shear_potential_vorticity_surface'
%     'U-Component_Storm_Motion_height_above_ground_layer'
%     'V-Component_Storm_Motion_height_above_ground_layer'
%     'Pressure_maximum_wind'
%     'Pressure_tropopause'
%     'Pressure_msl'
%     'Pressure_potential_vorticity_surface'
%     'Pressure_reduced_to_MSL_msl'
%     'Geopotential_height_zeroDegC_isotherm'
%     'Geopotential_height_maximum_wind'
%     'Geopotential_height_tropopause'
%     'Geopotential_height_isobaric'
%     'Geopotential_height_potential_vorticity_surface'
%     'Geopotential_height_highest_tropospheric_freezing'
%     'Geopotential_height_anomaly_isobaric'
%     '5-Wave_Geopotential_Height_isobaric'
%     '5-Wave_Geopotential_Height_Anomaly_isobaric'
%     'Cloud_water_entire_atmosphere'
%     'Parcel_lifted_index_to_500_hPa_pressure_difference_layer'
%     'Convective_available_potential_energy_surface'
%     'Convective_available_potential_energy_pressure_difference_layer'
%     'Convective_inhibition_surface'
%     'Convective_inhibition_pressure_difference_layer'
%     'Storm_relative_helicity_height_above_ground_layer'
%     'Surface_Lifted_Index_surface'
%     'Best_4_layer_Lifted_Index_surface'
%     'Total_ozone_entire_atmosphere'
%     'Ozone_Mixing_Ratio_isobaric'
%     'lat'
%     'lon'
%     'isobaric'
%     'isobaric1'
%     'pressure_difference_layer'
%     'height_above_ground'
%     'pressure_difference_layer1'
%     'potential_vorticity_surface'
%     'isobaric2'
%     'sigma'
%     'sigma_layer'
%     'height_above_ground_layer'
%     'altitude_above_msl'
%     'height_above_ground_layer1'
%     'isobaric3'
%     'time'
%     'time1'