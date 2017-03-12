function [v_grid, CIA_spectrum] = F_CIA2Spec(gfit_cia,Temperature,Pressure,vStart,vEnd,MR,v_grid)
% This function transform the GFIT psedo line list data to real spectrum
% using HITRAN like voigt profile

% Written by Kang Sun on 2016/7/20

% speed of light in SI unit
c = 2.99792458e8;
% Planck constant in SI unit
h = 6.62607004e-34;
% Bolzmann constant in SI unit
kB = 1.38064852e-23;
% Avogadro's Constant in SI unit
Na = 6.02214e23;
% second radiation constant, 1.4388 cm K, note it is cm
c2 = h*c/kB*100;

% HITRAN reference temperature/pressure, 296 K, 1 atm
T0 = 296; P0 = 1013.25;

Filter = gfit_cia.moleculeNumber == 7 & ...
    gfit_cia.transitionWavenumber <= vEnd+20 &...
    gfit_cia.transitionWavenumber >= vStart-20;

v0 = gfit_cia.transitionWavenumber(Filter);
S0 = gfit_cia.lineIntensity(Filter);
GammaP0 = gfit_cia.airBroadenedWidth(Filter);
GammaPs0 = gfit_cia.selfBroadenedWidth(Filter);
E = gfit_cia.lowerStateEnergy(Filter);
n = gfit_cia.temperatureDependence(Filter);
delta = gfit_cia.pressureShift(Filter);

%         % number densities in molecules/cm3
%         N_vec = MR*P_layer(ilayer)*100/kB/Temperature/1e6;
% Doppler HWHM at this layer
GammaD = (v0/c).*sqrt(2*kB*Na*Temperature*log(2)./0.032);
% line strength at this layer
S = S0*Q(T0,'O2')...
    /Q(Temperature,'O2')...
    .*exp(-c2*E/Temperature)./exp(-c2*E/T0)...
    .*(1-exp(-c2*v0/Temperature))./(1-exp(-c2*v0/T0));
% Lorentzian HWHM at this layer
GammaP = (GammaP0.*(1-MR)+GammaPs0.*MR)...
    *(Pressure/P0)...
    .*((T0/Temperature).^n);
% line position currected by air-broadened pressure shift
v0 = v0+delta*Pressure/P0;
% the core of line-by-line calculation
% for pseudo line list, v grid is identical to line centers
if ~exist('v_grid','var')
v_grid = v0;
end
CIA_spectrum = zeros(size(v_grid));
for iline = 1:length(v0)
    vFilter = v_grid > v0(iline)-50 ...
        & v_grid < v0(iline)+50;
    lineprofile = voigt_fn_fast(v_grid(vFilter),...
        v0(iline),GammaD(iline),GammaP(iline));
    temp = zeros(size(v_grid));
    temp(vFilter) = S(iline).*lineprofile;
    CIA_spectrum = CIA_spectrum+temp;
end

