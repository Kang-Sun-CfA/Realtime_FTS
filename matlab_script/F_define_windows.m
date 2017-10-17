function [vStart, vEnd, wStart, wEnd, mol_for_fit,mol_for_spec] = ...
    F_define_windows(target_gas,windowID)
% This function outputs parameters for different fitting windows
% mol_for_fit is name of mol included in the fitting,
% should include isotope like 'HDO'

% mol_for_spec is cell of strings of mol for spectroscopy calculation, no
% need for isotope. 'H2O' includes 'HDO' already.

% vStart, vEnd are the wave number boundaries of high res tau calculation

% wStart, wEnd are the wave number boudnaries of obs window, a little
% smaller than vStart and vEnd

% Written by Kang Sun on 2016/07/18
% Updated on 2016/08/02 to remove separate cntm, remove comma separated mol
% list
% Updated on 2017/10/07 to add O2 1.06 micron window for CIA
if ~exist('windowID','var')
    windowID = 1;
end
switch target_gas
    case 'CO2'
        switch windowID
            case 1
                wStart = 6220-40; wEnd = 6220+40;
                vStart = wStart-20; vEnd = wEnd+20;
                mol_for_fit = {'CO2','H2O','CH4'}; %HDO
                mol_for_spec = {'CO2','H2O','CH4'};
            case 2
                wStart = 6339.5-85/2; wEnd = 6339.5+85/2;
                vStart = wStart-20; vEnd = wEnd+20;
                mol_for_fit = {'CO2','H2O','HDO'};%HDO
                mol_for_spec = {'CO2','H2O'};
        end
    case 'CH4'
        switch windowID
            case 1
                wStart = 5880; wEnd = 5996;
                vStart = wStart-20; vEnd = wEnd+20;
                mol_for_fit = {'CH4','H2O'};
                mol_for_spec = {'CH4','H2O','CO2','N2O'};
            case 2
                wStart = 6007; wEnd = 6145;
                vStart = wStart-20; vEnd = wEnd+20;
                mol_for_fit = {'CH4','H2O','CO2'};
                mol_for_spec = {'CH4','H2O','CO2'};
        end
    case 'O2'
        switch windowID
            case 1
                vStart = 7700; vEnd = 8050;
                wStart = 7765; wEnd = 8005;
                mol_for_fit = {'O2','H2O'};
                mol_for_spec = {'H2O','O2'};
            case 0
                vStart = 9100; vEnd = 9700;
                wStart = 9120; wEnd = 9680;
                mol_for_fit = {'O2','H2O'};
                mol_for_spec = {'H2O','O2'};
        end
    case 'H2O'
        switch windowID
            case 1
                wStart = 6076-3.85/2;wEnd = 6076+3.85/2;
                vStart = wStart-10; vEnd = wEnd+10;
                mol_for_fit = {'H2O','CH4','HDO','CO2'};
                mol_for_spec = {'CH4','H2O','CO2'};
        end
end
