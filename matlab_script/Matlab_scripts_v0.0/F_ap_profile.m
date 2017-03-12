function profile = F_ap_profile(taumol,profileinput)
% function to generate a priori profiles based on forecast data or GFIT vmr
% CO2 is adjusted to 390 ppm, CH4 is adjusted to 1.9 ppm
C = profileinput.C;
Z_layer = profileinput.Z_layer;
cfsH2O = profileinput.cfsH2O;

switch taumol
    case 'H2O'
        profile = cfsH2O;
    case 'CO2'
        profile = interp1(C(:,1)*1000,C(:,3),Z_layer);
        proflie = profile/profile(1)*390e-6;
    case 'O3'
        profile = interp1(C(:,1)*1000,C(:,4),Z_layer);
    case 'N2O'
        profile = interp1(C(:,1)*1000,C(:,5),Z_layer);
    case 'CO'
        profile = interp1(C(:,1)*1000,C(:,6),Z_layer);
    case 'CH4'
        profile = interp1(C(:,1)*1000,C(:,7),Z_layer);
        proflie = profile/profile(1)*1.9e-6;
    case 'O2'
        profile = interp1(C(:,1)*1000,C(:,8),Z_layer);
    otherwise
        profile = nan(size(Z_layer));
end