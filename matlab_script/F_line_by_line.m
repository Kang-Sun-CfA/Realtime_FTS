function output_LBL = F_line_by_line(input_LBL)
% this function calcuates the layered optical depth for a single molecule.
% Updated by Kang Sun on 2016/7/16
% Updated by Kang Sun on 2016/7/30 to add savepath

% Major revision by Kang Sun on 2016/8/2 to use common grid and output as
% layered matrices, not cells. All layers are defined at the common grid
% Basically, F_line_by_line and F_generate_tau_mat are combined, except:
% Line by line calc is made at 3 per HWHM, then convolved to 3 x
% common_grid interval, then resampled to commong_grid

% all above were from the older version F_line_by_line_FTS.m

% Resaved as F_line_by_line.m by Kang Sun on 2017/01/24 to simplify inputs
% and add options of different CIA

% Revised by Kang Sun on 2017/02/08 to add more CIA options

taumol = input_LBL.taumol;
profile = input_LBL.profile;
P_level = input_LBL.P_level;
P_layer = input_LBL.P_layer;
T_layer = input_LBL.T_layer;
Z_level = input_LBL.Z_level;
Z_layer = input_LBL.Z_layer;
common_grid = input_LBL.common_grid;
vStart = min(common_grid);
vEnd = max(common_grid);
molparampath = input_LBL.molparampath;
datadir = input_LBL.datadir;
if isfield(input_LBL,'which_CIA')
    which_CIA = input_LBL.which_CIA;
else
    which_CIA = 'gfit';
end

if ~iscell(which_CIA)
    which_CIA = {which_CIA};
end
hitranpath = [datadir,'FTS_5000to12000','.mat'];

% load O2 CIA at 1.27 um band
if strcmpi(taumol,'O2') && mean(common_grid) < 8500
    if sum(ismember(which_CIA,'gfit'))
        fciapath = [datadir,'gfit_fcia.mat'];
        sciapath = [datadir,'gfit_scia.mat'];
        
        temp = load(fciapath,'gfit_fcia');
        gfit_fcia = temp.gfit_fcia;
        temp = load(sciapath,'gfit_scia');
        gfit_scia = temp.gfit_scia;
        
        v_grid_gfit = vStart:0.02:vEnd;
        vlow_gfit = vStart:1:vEnd;
        [~,~,bin] = histcounts(v_grid_gfit,vlow_gfit);
    end
    if sum(ismember(which_CIA,'sao'))
        %     % these are older files given by Tijs on 2017/01/13
        %     fid = fopen([datadir,'av0.table']);
        %     C_sao = cell2mat(textscan(fid,'%f%f%f%f','headerlines',1,'delimiter',' ',...
        %         'multipledelimsasone',1));
        %     vlow_sao = C_sao(:,1);
        %     [xgrid,ygrid] = meshgrid(vlow_sao,[298 273 253]);
        
        % these are newer files given by Tijs on 2017/01/27
        fid = fopen([datadir,'ani_127.table']);
        C_sao = cell2mat(textscan(fid,'%f%f%f%f%f%f%f','headerlines',1,'delimiter',' ',...
            'multipledelimsasone',1));
        vlow_sao = C_sao(:,1);
        cia_temp_vec = [298 273 253 228 203 178];
        [xgrid,ygrid] = meshgrid(vlow_sao,cia_temp_vec);
    end
    if sum(ismember(which_CIA,'mate'))
        ciafname1 = [datadir,'O2-Air_Mate.cia'];
        % ciafname1 = [datadir,'O2-Air_SN_126.cia'];
        C = cell(2,1);
        headers = cell(2,1);
        lineN = nan(2,1);
        ciaT = nan(2,1);
        fid = fopen(ciafname1);
        count = 0;
        while 1
            count = count+1;
            tline = fgetl(fid);
            if ~ischar(tline)
                break;
            end
            headers{count} = tline;
            temp = textscan(tline,'%s%f%f%f%f%f%f','Delimiter',' ','MultipleDelimsAsOne',1);
            lineN(count) = temp{4};
            ciaT(count) = temp{5};
            C(count) = textscan(fid,'%f%f','Delimiter',' ','MultipleDelimsAsOne',1,...
                'collectoutput',1);
        end
        fclose(fid);
        ciaT_mate = ciaT;
        vlow_mate = C{1}(:,1);
        C_mate = repmat(C{1}(:,2),[1,length(C)]);
        for imate = 2:length(C)
            C_mate(:,imate) = interp1(C{imate}(:,1),C{imate}(:,2),vlow_mate,'linear','extrap');
        end
        [xgrid_mate,ygrid_mate] = meshgrid(vlow_mate,[253 273 296]);
    end
end

% load O2 CIA at 1.06 um band
if strcmpi(taumol,'O2') && mean(common_grid) >= 8500
    if sum(ismember(which_CIA,'gfit'))
        fciapath = [datadir,'gfit_fcia.mat'];
        sciapath = [datadir,'gfit_scia.mat'];
        
        temp = load(fciapath,'gfit_fcia');
        gfit_fcia = temp.gfit_fcia;
        temp = load(sciapath,'gfit_scia');
        gfit_scia = temp.gfit_scia;
        
        v_grid_gfit = vStart:0.02:vEnd;
        vlow_gfit = vStart:1:vEnd;
        [~,~,bin] = histcounts(v_grid_gfit,vlow_gfit);
    end
    if sum(ismember(which_CIA,'sao'))
        fid = fopen([datadir,'iso_106.table']);
        C_sao = cell2mat(textscan(fid,'%f%f%f%f%f%f%f','headerlines',1,'delimiter',' ',...
            'multipledelimsasone',1));
        vlow_sao = C_sao(:,1);
        cia_temp_vec = [298 273 253 228 203 178];
        [xgrid,ygrid] = meshgrid(vlow_sao,cia_temp_vec);
    end
    C_mate = [];xgrid_mate = [];ygrid_mate = [];vlow_mate = [];
end
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

mol_list = {'H2O','CO2','O3','N2O','CO','CH4','O2','HF'};
% all 1 if only look at the major isotope
nisotope = ones(length(mol_list),1);
% consider HDO
nisotope(1) = 4;
% read molar weight, mol number and abundances from molparam.txt
M_list = loadmolparam(molparampath,mol_list,nisotope);

dZ = abs(diff(Z_level));
nlevel = length(P_level);
nlayer = nlevel-1;

% if, O2, add CIA, use GFIT pseudo line list
if strcmpi(taumol,'O2')
    if sum(ismember(which_CIA,'gfit'))
        dtau_cntm_gfit = zeros(length(common_grid),nlayer,'single');
    end
    if sum(ismember(which_CIA,'sao'))
        dtau_cntm_sao = zeros(length(common_grid),nlayer,'single');
    end
    if sum(ismember(which_CIA,'mate'))
        dtau_cntm_mate = zeros(length(common_grid),nlayer,'single');
    end
    parfor ilayer = 1:nlayer
        MR_O2 = profile(ilayer); % has to be profile of O2, 0.2095
        MR_N2 = 0.78;
        % number density of air in this layer in molec/cm3
        N_ilayer = P_layer(ilayer)*100/1e6/kB/T_layer(ilayer);
        if sum(ismember(which_CIA,'gfit'))
            [~, fcia_temp] = ...
                F_CIA2Spec(gfit_fcia,T_layer(ilayer),P_layer(ilayer),vStart,vEnd,1,v_grid_gfit);
            [~, scia_temp] = ...
                F_CIA2Spec(gfit_scia,T_layer(ilayer),P_layer(ilayer),vStart,vEnd,1,v_grid_gfit);
            fcia_low = accumarray(bin(:),fcia_temp,[],@mean);
            scia_low = accumarray(bin(:),scia_temp,[],@mean);
            acia_low = MR_O2*scia_low + MR_N2*fcia_low;
            dtau_cntm_gfit(:,ilayer) = N_ilayer^2*MR_O2*100*dZ(ilayer)*...
                interp1(0.5*(vlow_gfit(1:end-1)+vlow_gfit(2:end)),acia_low,common_grid);
        end
        if sum(ismember(which_CIA,'sao'))
            if T_layer(ilayer) >= max(cia_temp_vec)
                acia_low = C_sao(:,2);
            elseif T_layer(ilayer) <= min(cia_temp_vec)
                acia_low = C_sao(:,end);
            else
                acia_low = interp2(xgrid,ygrid,C_sao(:,2:end)',vlow_sao,T_layer(ilayer));
            end
            dtau_cntm_sao(:,ilayer) = N_ilayer^2*MR_O2*100*dZ(ilayer)*...
                interp1(vlow_sao,acia_low,common_grid);
        end
        if sum(ismember(which_CIA,'mate'))
            if T_layer(ilayer) >= 296
                acia_low = C_mate(:,3);
            elseif T_layer(ilayer) <= 253
                acia_low = C_mate(:,1);
            else
                acia_low = interp2(xgrid_mate,ygrid_mate,C_mate',vlow_mate,T_layer(ilayer));
            end
            dtau_cntm_mate(:,ilayer) = N_ilayer^2*MR_O2*100*dZ(ilayer)*...
                interp1(vlow_mate,acia_low,common_grid);
        end
    end
    if sum(ismember(which_CIA,'sao'))
        output_LBL.dtau_cntm_sao = dtau_cntm_sao;
    end
    if sum(ismember(which_CIA,'gfit'))
        output_LBL.dtau_cntm_gfit = dtau_cntm_gfit;
    end
    if sum(ismember(which_CIA,'mate'))
        output_LBL.dtau_cntm_mate = dtau_cntm_mate;
    end
end

%% construct the wavenumber grid
nsample = 3; % how many samples per HWHM?
rep_HWHM = nan(nlayer,3);
v_grid_cell = cell(nlayer,1);
for ilayer = 1:nlayer
    % imagine a molecule with the lorentizian width of methane and gaussian
    % width of O3, assume n = 0.5 (whatever). This is because the line-
    % intensity-weighted pressure broadening of methane is the smallest
    % ammong major GHGs and the gaussian width of O3 is the smallest.
    % Assume the lorentzian width is constant at the weighted mean, the
    % gaussian width is constant at the value of 1000 cm-1
    rep_HWHML = 0.03*(P_layer(ilayer)/P0).*((T_layer(ilayer)/T0).^0.5);
    rep_HWHMD = mean([vStart vEnd])/c*sqrt(2*kB*Na*T_layer(ilayer)*log(2)/0.044);
    % calculate voigt width using lorentzian and gaussian width, using
    % equations 5-6 from Olivero, J. J., and R. L. Longbothum. "Empirical
    % fits to the Voigt line width: A brief review."
    d = (rep_HWHML-rep_HWHMD)/(rep_HWHML+rep_HWHMD);
    Beta = 0.023665*exp(0.6*d)+0.00418*exp(-1.9*d);
    rep_HWHMV = (1-0.18121*(1-d^2)-Beta*sin(pi*d))*(rep_HWHML+rep_HWHMD);
    rep_HWHM(ilayer,1) = rep_HWHML;rep_HWHM(ilayer,2) = rep_HWHMD;
    rep_HWHM(ilayer,3) = rep_HWHMV;
    % quite weird, the start:step:end construction method gives very
    % inaccurate diff(v_grid) values
    % sorry single is not precise enough for such exotic v grid, need
    % double precision
    %     if rep_HWHMV/nsample >= median(diff(common_grid))
    v_grid_cell{ilayer} = ...
        linspace(double(vStart),double(vEnd),...
        double(round((vEnd-vStart)/(rep_HWHMV/nsample))));
    %     else
    %         v_grid_cell{ilayer} = common_grid;
    %     end
    
end
% cutoff = rep_HWHM2(:,3)*1000;
cutoff = ones(size(Z_layer))*25;
% loglog(cutoff,P_layer);set(gca,'yscale','log','ydir','rev')
% semilogx(rep_HWHM/nsample,Z_layer)
output_LBL.rep_HWHM = rep_HWHM;
%% load hitran file
Swithlines = load(hitranpath,'lines');
%% calculate layer optical depth, slowest part
% tic
if isempty(taumol) % if no input, calculate N2O, because it's fast
    includemol_list = {'N2O'};
else
    includemol_list = {taumol};
end

% ntime = nan(nlayer,1);
localmollist = mol_list(:);

idx = strcmp(localmollist, taumol);
howmanyiso = nisotope(idx);
dtau_layer = nan(howmanyiso,length(common_grid),nlayer,'single');

parfor ilayer = 1:nlayer
    %     tic
    % degrade tau of each layer to this FWHM:
    dgrd_fwhm = 2.5*median(diff(common_grid));
    v_grid = v_grid_cell{ilayer};
    localMlist = M_list;
    lines = Swithlines.lines;
    localnisotope = nisotope;
    for imol = 1:length(includemol_list)
        
        idx = find(strcmp(localmollist, includemol_list{imol}));
        
        molnumber = localMlist(idx).molnumber;
        mixingratio = profile(ilayer);
        howmanyiso = localnisotope(idx);
        % initialize tau for each species
        tau_grid = zeros(length(v_grid),howmanyiso,'single');
        for iiso = 1:howmanyiso
            MW = localMlist(idx).isotopologues(iiso).molarmass;
            abundance = localMlist(idx).isotopologues(iiso).abundance;
            Filter = lines.moleculeNumber == molnumber & ...
                lines.isotopologueNumber == iiso & ...
                lines.transitionWavenumber >= vStart-12.5 &...
                lines.transitionWavenumber <= vEnd+12.5;
            isoN = lines.isotopologueNumber(Filter);
            v0 = lines.transitionWavenumber(Filter);
            S0 = lines.lineIntensity(Filter);
            GammaP0 = lines.airBroadenedWidth(Filter);
            GammaPs0 = lines.selfBroadenedWidth(Filter);
            E = lines.lowerStateEnergy(Filter);
            n = lines.temperatureDependence(Filter);
            delta = lines.pressureShift(Filter);
            MR = mixingratio*abundance; % no consider water vapor dillution
            
            % number densities in molecules/cm3
            N = MR*P_layer(ilayer)*100/kB/T_layer(ilayer)/1e6;
            % Doppler HWHM at this layer
            GammaD = (v0/c).*sqrt(2*kB*Na*T_layer(ilayer)*log(2)./MW);
            % line strength at this layer
            S = S0*Q(T0,includemol_list{imol},iiso)...
                /Q(T_layer(ilayer),includemol_list{imol},iiso)...
                .*exp(-c2*E/T_layer(ilayer))./exp(-c2*E/T0)...
                .*(1-exp(-c2*v0/T_layer(ilayer)))./(1-exp(-c2*v0/T0));
            % Lorentzian HWHM at this layer
            GammaP = (GammaP0.*(1-MR)+GammaPs0.*MR)...
                *(P_layer(ilayer)/P0)...
                .*((T0/T_layer(ilayer)).^n);
            % line position currected by air-broadened pressure shift
            v0 = v0+delta*P_layer(ilayer)/P0;
            % the core of line-by-line calculation
            localcutoff = cutoff(ilayer);
            localdZ = dZ(ilayer);
            for iline = 1:length(v0)
                vFilter = v_grid > v0(iline)-localcutoff...
                    & v_grid < v0(iline)+localcutoff;
                lineprofile = voigt_fn_fast(v_grid(vFilter),...
                    v0(iline),GammaD(iline),GammaP(iline));
                if molnumber == 1 && iiso == 1 % remove baseterm for water vapor
                    lineprofile = lineprofile-min(lineprofile);
                end
                localtau = zeros(length(v_grid),1,'single');
                localtau(vFilter) = N*100*localdZ*S(iline).*lineprofile;
                tau_grid(:,iiso) = tau_grid(:,iiso)+localtau;
                %             tau_grid(vFilter) = tau_grid(vFilter)+N_vec(iline)...
                %                 .*100*localdZ*S(iline).*lineprofile;
            end
        end
        tau_common_grid = zeros(length(common_grid),howmanyiso,'single');
        
        for iiso = 1:howmanyiso
            tau_common_grid(:,iiso) =...
                F_conv_interp(v_grid,tau_grid(:,iiso),dgrd_fwhm,common_grid)
        end
        dtau_layer(:,:,ilayer) = tau_common_grid';
    end
    %     ntime(ilayer) = toc;
    %     disp(['layer ',num2str(ilayer),' is finished at ',datestr(now),...
    %         ', took ',num2str(ntime(ilayer)),' s.'])
    
    
end
% kind of creepy. take only HDO out. isotope 2 and 3 are useless anyway
if strcmp(taumol,'H2O')
    dtau_layer = dtau_layer([1 4],:,:);
end
output_LBL.dtau_layer = dtau_layer;
output_LBL.taumol_info = M_list(strcmp(mol_list(:), includemol_list{1}));
% output_LBL.common_grid = common_grid;


function M_list = loadmolparam(filename,mol_list,nisotope)
if nargin < 3
    nisotope = ones(length(mol_list),1);
end
% Read the molecular parameters file into a character
% array.
text = fileread(filename);

% Get the list of molecules.
molecules = regexp(text, '^\s*([\w\+\-\d]+)\s+\(\d+\)\s*\n', 'lineanchors', 'tokens');

% Split the file into sections: one for each molecule.
splitStrings = regexp(text, '^\s*[\w\+\-\d]+\s+\(\d+\)\s*\n', 'split', 'lineanchors');
molDataText = splitStrings(2:end); % Throw away the first line and ensure there's a newline on the end.
numberOfMolecules = numel(molDataText);

M_list = struct('formula',{},'molnumber',{},'isotopologues',{});
countm = 1;
for m = 1:numberOfMolecules
    % Scan the text for the isotopologue information and store it
    % in the matrix I.
    isoData = textscan(molDataText{m}, '%f %f %f %f %f');
    %     numberOfIsotopologues = numel(isoData{1});
    % The following line is here because of typing problems.
    if countm <= length(mol_list) && strcmp(mol_list{countm},char(molecules{m}))
        M_list(countm).formula = char(molecules{m});
        M_list(countm).molnumber = m;
        % Copy the scanned data into a structure.
        for n = 1:nisotope(countm)
            M_list(countm).isotopologues(n).number = isoData{1}(n);
            M_list(countm).isotopologues(n).abundance = isoData{2}(n);
            M_list(countm).isotopologues(n).Q = isoData{3}(n);
            M_list(countm).isotopologues(n).gj = isoData{4}(n);
            M_list(countm).isotopologues(n).molarmass = isoData{5}(n)/1000;
        end
        countm = countm+1;
    end    
end