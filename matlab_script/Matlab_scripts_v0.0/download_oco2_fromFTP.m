%% download data from OCO-2 ftp
% clear
clc
cd('/data/tempo1/Shared/kangsun/oco-2')
addpath('/home/kangsun/Solar')
addpath(genpath('/home/kangsun/matlab functions'))
f = ftp('oco2.gesdisc.eosdis.nasa.gov:21');

% oYear= 2015;
% oDOYvec = 335:366;%-116:335;

oDOYvec = 367:453;%-116:335;
for ioDOY = 1:length(oDOYvec)
    if oDOYvec(ioDOY) <= 0
        oYear = 2014;
    elseif oDOYvec(ioDOY) > 365
        oYear = 2016;
    else
        oYear = 2015;
    end
DOY = mod(oDOYvec(ioDOY),365);
if DOY == 0
DOY = 365;
end
%Dirname = ['data/s4pa/OCO2_DATA/OCO2_L2_Diagnostic.7r/',num2str(oYear),'/',num2str(oDOYvec(ioDOY),'%.3d'),'/'];
    Dirname = ['data/s4pa/OCO2_DATA/OCO2_L1B_Calibration.7r/',num2str(oYear),'/',num2str(DOY,'%.3d'),'/'];
    details = dir(f,Dirname);
    if ~isempty(details)
        for ifile = 1:length(details)
        %ifile = 0; Criterion = 0;
        %while Criterion == 0
         %   ifile = ifile+1;
            %     if ~isempty(regexp(details(ifile).name,'L2DiaGL', 'once')) && isempty(regexp(details(ifile).name,'xml', 'once'))
            if ~isempty(regexp(details(ifile).name,'L1bClSS', 'once')) && isempty(regexp(details(ifile).name,'xml', 'once'))
                mget(f,[Dirname,details(ifile).name]);
                
          %      Criterion = 1;
            end
        end
    end
end
close(f)
%% find SDD
clc
oDOYvec = 367:453;%-116:335;
for ioDOY = 1:length(oDOYvec)
            if oDOYvec(ioDOY) <= 0
                oYear = 2014;
            elseif oDOYvec(ioDOY) > 365
                oYear = 2016;
            else
                oYear = 2015;
            end
            DOY = mod(oDOYvec(ioDOY),365);
            if DOY == 0
                DOY = 365;
            end
            Dirname = ['/data/tempo1/Shared/kangsun/oco-2/data/s4pa/OCO2_DATA/',...
                'OCO2_L1B_Calibration.7r/',num2str(oYear),'/',...
                num2str(DOY,'%.3d'),'/'];
            details = dir(Dirname);
            if ~isempty(details)
                %         nssfile = 0;
                for ifile = 1:length(details)
                    
                    if details(ifile).bytes > 100000000
                        disp(details(ifile).name)
                        disp(Dirname)
                    end
                end
            end
end
%% download scientific L1B data
% clear
clc
cd('/data/tempo1/Shared/kangsun/oco-2')
addpath('/home/kangsun/Solar')
addpath(genpath('/home/kangsun/matlab functions'))
f = ftp('oco2.gesdisc.eosdis.nasa.gov:21');

% oYear= 2015;
% oDOYvec = 335:366;%-116:335;

oDOYvec = 316;%-116:335;
for ioDOY = 1:length(oDOYvec)
    if oDOYvec(ioDOY) <= 0
        oYear = 2014;
    elseif oDOYvec(ioDOY) > 365
        oYear = 2016;
    else
        oYear = 2015;
    end
DOY = mod(oDOYvec(ioDOY),365);
if DOY == 0
DOY = 365;
end
%Dirname = ['data/s4pa/OCO2_DATA/OCO2_L2_Diagnostic.7r/',num2str(oYear),'/',num2str(oDOYvec(ioDOY),'%.3d'),'/'];
    Dirname = ['data/s4pa/OCO2_DATA/OCO2_L1B_Science.7r/',num2str(oYear),'/',num2str(DOY,'%.3d'),'/'];
    details = dir(f,Dirname);
    if ~isempty(details)
        for ifile = 1:length(details)
        %ifile = 0; Criterion = 0;
        %while Criterion == 0
         %   ifile = ifile+1;
            %     if ~isempty(regexp(details(ifile).name,'L2DiaGL', 'once')) && isempty(regexp(details(ifile).name,'xml', 'once'))
            if ~isempty(regexp(details(ifile).name,'L1bSc', 'once')) && isempty(regexp(details(ifile).name,'xml', 'once'))
                mget(f,[Dirname,details(ifile).name]);
                
          %      Criterion = 1;
            end
        end
    end
end
close(f)
                        