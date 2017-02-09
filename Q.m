function Q = Q(T,mol,iiso)

% Reference: Gamache et al.,(2000)
% Total internal partition sums for molecules in the terrestrial atmosphere

if nargin < 2
    mol = 'default';iiso = 1;
elseif nargin < 3
    iiso = 1;
end

switch mol
    case 'NO'
        p = [0 0 0 1];
    case 'N2O'
        p = [0.46310e-4 -0.76213e-2 0.14979e2 0.24892e2];
    
    case 'CH4'
        p = [0.15117e-5 0.26831e-2 0.11557e1 -0.26479e2];
         
    case 'CO'
        p = [0.14896e-7 -0.74669e-5 0.3629 0.27758];
        
    case 'C2H2'
        p = [0.84612e-5 -0.25946e-2 0.14484e1 -0.83088e1];
        
    case 'NH3'
        p = [0.18416e-5 0.94575e-2 0.30915e1 -0.62293e2];
        
    case 'C2H4'
        p = [0 0 0 1];
        
    case 'H2O'
        switch iiso
            case 1
                p = [0.48938e-6 0.12536e-2 0.27678 -0.44405e1];
            case 2
                p = [0.52046e-6 0.12802e-2 0.27647 -0.43624e1];
            case 3
                p = [0.31668e-5 0.76905e-2 0.16458e1 -0.25767e2];
            case 4
                p = [0.21530e-5 0.61246e-2 0.13793e1 -0.23916e2];
            otherwise
                p = [0.48938e-6 0.12536e-2 0.27678 -0.44405e1];
        end
        
    case 'CO2'
        p = [0.25974e-5 -0.69259e-3 0.94899 -0.1317e1];
        
    case 'O3'
        p = [0.26669e-4 0.10396e-1 0.69047e1 -0.16443e3];
        
    case 'SO2'
        p = [0.52334e-4 0.22164e-1 0.11101e2 -0.24056e3];
        
    case 'O2'
        p = [0.13073e-6 -0.64870e-4 0.73534 0.35923];
    
    case 'HF'
        p = [0.46889e-8 0.59154e-5 0.13350 0.15486e1];
        
    case 'default'
        p = [0 0 0 1];
    otherwise
        p = [0 0 0 1];
        
end

Q = polyval(p,T);