% calculate the Voigt profile using the complex error function

function lsSA2 = voigt_fn_fast(nu,nu0,gamma_D,gamma_L)

% nu is wavenumber array [cm^-1]
% nu0 is line center [cm^-1]
% gamma_D is Doppler linewidth [cm^-1]
% gamma_L is Lorentz linewidth [cm^-1]

% convert to dimensionless units
x = sqrt(log(2)).*(nu-nu0)./(gamma_D);
y = sqrt(log(2)).*(gamma_L/gamma_D);

% call complexErrorFunction
lsSA2 =(y/sqrt(pi)/gamma_L)*real(voigtf(x,y,2));

if size(lsSA2) ~= size(nu)
    lsSA2 = lsSA2';
end
end

function VF = voigtf(x,y,opt)

% This function file is a subroutine for computation of the Voigt function.
% The input parameter y is used by absolute value. The parameter opt is 
% either 1 for more accurate or 2 for more rapid computation. 
%
% NOTE: This program completely covers the domain 0 < x < 40,000 and 
% 10^-4 < y < 10^2 required for applications using HITRAN molecular 
% spectroscopic database. However, it may be implemented only to cover the 
% smaller domain 0 <= x <= 15 and 10^-6 <= y <= 15 that is the most 
% difficult for rapid and accurate computation.
%
% The code is written by Sanjar M. Abrarov and Brendan M. Quine, York
% University, Canada, March 2015.

if nargin == 2
    opt = 1; % assign the defaul value opt = 1
end

if opt ~= 1 && opt ~=2
    disp(['opt = ',num2str(opt),' cannot be assigned. Use either 1 or 2.'])
    return
end

% *************************************************************************
% Define array of coefficients as coeff = [alpha;beta;gamma]'
% *************************************************************************
if opt == 1
    
    coeff = [
    1.608290174437121e-001 3.855314219175531e-002  1.366578214428949e+000
    6.885967427017463e-001 3.469782797257978e-001 -5.742919588559361e-002
    2.651151642675390e-001 9.638285547938826e-001 -5.709602545656873e-001
   -2.050008245317253e-001 1.889103967396010e+000 -2.011075414803758e-001
   -1.274551644219086e-001 3.122804517532180e+000  1.069871368716704e-002
   -1.134971805306579e-002 4.664930205202391e+000  1.468639542320982e-002
    4.201921570328543e-003 6.515481030406647e+000  1.816268776500938e-003
    8.084740485193432e-004 8.674456993144942e+000 -6.875907999947567e-005
    1.946391440605860e-005 1.114185809341728e+001 -2.327910355924500e-005
   -4.132639863292073e-006 1.391768433122366e+001 -1.004011418729134e-006
   -2.656262492217795e-007 1.700193570656409e+001  2.304990232059197e-008
   -1.524188131553777e-009 2.039461221943855e+001  2.275276345355270e-009
    2.239681784892829e-010 2.409571386984707e+001  3.383885053101652e-011
    4.939143128687883e-012 2.810524065778962e+001 -4.398940326332977e-013
    4.692078138494072e-015 3.242319258326621e+001 -1.405511706545786e-014
   -2.512454984032184e-016 3.704956964627684e+001 -3.954682293307548e-016
    ]; 
    mMax = 16; % 16 summation terms

elseif opt == 2
        
    coeff = [
    2.307372754308023e-001 4.989787261063716e-002  1.464495070025765e+000
    7.760531995854886e-001 4.490808534957343e-001 -3.230894193031240e-001
    4.235506885098250e-002 1.247446815265929e+000 -5.397724160374686e-001
   -2.340509255269456e-001 2.444995757921221e+000 -6.547649406082363e-002
   -4.557204758971222e-002 4.041727681461610e+000  2.411056013969393e-002
    5.043797125559205e-003 6.037642585887094e+000  4.001198804719684e-003
    1.180179737805654e-003 8.432740471197681e+000 -5.387428751666454e-005
    1.754770213650354e-005 1.122702133739336e+001 -2.451992671326258e-005
   -3.325020499631893e-006 1.442048518447414e+001 -5.400164289522879e-007
   -9.375402319079375e-008 1.801313201244001e+001  1.771556420016014e-008
    8.034651067438904e-010 2.200496182129099e+001  4.940360170163906e-010
    3.355455275373310e-011 2.639597461102705e+001  5.674096644030151e-014
    ];
    mMax = 12; % 12 summation terms
end
% *************************************************************************

varsigma = 2.75; % define the shift constant
y = abs(y) + varsigma/2;

arr1 = y.^2 - x.^2; % define 1st repeating array
arr2 = x.^2 + y.^2; % define 2nd repeating array
arr3 = arr2.^2;  % define 3rd repeating array

    VF = 0; % initiate VF
    for m = 1:mMax
        VF = VF + (coeff(m,1)*(coeff(m,2) + arr1) + ...
            coeff(m,3)*y.*(coeff(m,2) + arr2))./(coeff(m,2)^2 + ...
            2*coeff(m,2)*arr1 + arr3);
    end
end
