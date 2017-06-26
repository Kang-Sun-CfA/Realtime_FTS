function window_list = F_absco2tau(inp)
% function used to generating spectroscopic fitting database offline, after
% the vertical profiles are available (e.g., from weather forecast)

% inputs: location of absco tables; vertical profiles saved in a correct
% format; common resolution (sampling interval really; the FWHM is 2.5 x
% common resolution)

% outputs: database for each species in each fitting window

% It has similar functionality to F_line_by_line.m, but interpolating ABSCO
% table instead of line-by-line voigt profile calculation.

% written by Kang Sun on 2017/06/20

% % inputs for testing
% clc;clear
% inp = [];
% inp.absco_dir = '/data/tempo1/Shared/kangsun/FTS_data/HITRAN/';
% inp.profile_fn = '~/FTS/toolbox/profiles_0000_0000.dat';
% inp.lowest_possible_Psurf = 980;
% inp.nsublayer = 3;
% inp.common_grid_resolution = 0.01;
% inp.surface_layer_P = 990;

% pressure of lowest level, below which is the adjustable surface layer
if ~isfield(inp,'lowest_possible_Psurf')
    inp.lowest_possible_Psurf = 980;
end

% representative pressure of the surface layer, should be close to the mean
% of inp.lowest_possible_Psurf and measured instantaneous surface pressure
if ~isfield(inp,'surface_layer_P')
    inp.surface_layer_P = 990;
end

% wavenumber interval to define the output
if ~isfield(inp,'common_grid_resolution')
    inp.common_grid_resolution = 0.05;
end

% how many sublayer to divide each layer into? OCO-2 used 10
if ~isfield(inp,'nsublayer')
    inp.nsublayer = 3;
end
varname = {'place_holder','Pressure','Temperature','Wavenumber'};

% Bolzmann constant in SI unit
kB = 1.38064852e-23;

dgrd_fwhm = 2.5*inp.common_grid_resolution;
%%
fid = fopen(inp.profile_fn);
C_0 = cell2mat(textscan(fid,'%f%f%f%f%f%f%f','delimiter',' ','multipledelimsasone',1,'headerlines',11));
fclose(fid);

% convert presgrid from Pa (or atm) to hPa if necessary
if max(C_0(:,1)) > 7.7e4 % I guess a good pressure profile extends below 770 hPa
    C_0(:,1) = C_0(:,1)/100;
end
if max(C_0(:,1)) < 2
    C_0(:,1) = C_0(:,1)*1.01325e3;
end

% make sure from low to high pressure, VERY IMPORTANT
[~,I] = sort(C_0(:,1));C_0 = C_0(I,:);

% trim the lowest part of profiles
C_low = nan(1,size(C_0,2));
C_low(1) = inp.lowest_possible_Psurf;
for i = 2:size(C_0,2)
    C_low(i) = interp1(C_0(:,1),C_0(:,i),inp.lowest_possible_Psurf);
end

C = cat(1,C_0(C_0(:,1)-inp.lowest_possible_Psurf < 0,:),C_low);

nsublayer = inp.nsublayer;
nlevel = size(C,1);
nlayer = nlevel-1;

%%%%%% Define retrieval windows
field1 = 'target_gas';
field2 = 'windowID';
field3 = 'common_grid';
field4 = 'tau_struct';
field5 = 'vRange';
field6 = 'wRange';
value1 = {'O2','CO2','CO2','CH4','CH4'};
value2 = {1,1,2,1,2};
window_list = struct(field1,value1,field2,value2,field3,[],field4,[],...
    field5,[],field6,[]);

for iwin = 1:length(window_list)
    target_gas = window_list(iwin).target_gas;
    windowID = window_list(iwin).windowID;
    [vStart, vEnd, wStart, wEnd, mol_for_fit] = ...
        F_define_windows(target_gas,windowID);
    common_grid = vStart:inp.common_grid_resolution:vEnd;
    window_list(iwin).common_grid = common_grid;
    window_list(iwin).vRange = [vStart,vEnd];
    window_list(iwin).wRange = [wStart,wEnd];
    tmp_struct = [];
    for imol = 1:length(mol_for_fit)
        absco_fn = [inp.absco_dir,target_gas,'_win',num2str(windowID),'_',...
            mol_for_fit{imol},'.h5'];
        if ~exist(absco_fn,'file')
            error(['No absco file ',target_gas,'_win',num2str(windowID),'_',...
                mol_for_fit{imol},'.h5!!!'])
        end
        switch mol_for_fit{imol}
            case 'O2'
                molN = 7;prof_index = 7;
            case 'CO2'
                molN = 2;prof_index = 5;
            case 'CH4'
                molN = 6;prof_index = 6;
            case {'H2O','HDO'}
                molN = 1;prof_index = 4;
        end
        varname{1} = ['Gas_',sprintf('%02d',molN),'_Absorption'];
        % read absco data
        absco = F_read_hdf(absco_fn,varname);
        inp_interp = [];
        inp_interp.tempgrid = absco.Temperature.data;
        inp_interp.presgrid = absco.Pressure.data;
        inp_interp.wavegrid = absco.Wavenumber.data;
        inp_interp.absco = absco.(varname{1}).data;
        inp_interp.Wq = absco.Wavenumber.data(absco.Wavenumber.data>=vStart & ...
            absco.Wavenumber.data<=vEnd);
        % optical depth of each layer, defined at common grid
        Tau_layer = zeros(nlayer,length(common_grid),'single');
        count_interp = 0;
        
        for ilevel = 1:nlevel-1
            C1 = C(ilevel,:);C2 = C(ilevel+1,:);
            dP = C2(1)-C1(1);
            dZ = -(C2(2)-C1(2));
            Psublayer = linspace(C1(1)+0.5*dP/nsublayer,C2(1)-0.5*dP/nsublayer,nsublayer);
            Psublevel = linspace(C1(1),C2(1),nsublayer+1);
            Zsublevel = interp1([C1(1) C2(1)],[C1(2) C2(2)],Psublevel);
            Csublayer = nan(nsublayer,size(C,2));
            Csublayer(:,1) = Psublayer(:);
            Csublayer(:,2) = -diff(Zsublevel);
            for iprof = 3:size(C,2);
                Csublayer(:,iprof) = interp1([C1(1) C2(1)],[C1(iprof) C2(iprof)],Psublayer);
            end
            
            % calculate optical depth of each layer
            dTau_layer = zeros(1,length(inp_interp.Wq),'single');
            for isublayer = 1:nsublayer
                eVMR = Csublayer(isublayer,prof_index);  % sublayer VMR
                % not sure about HDO
%                 if strcmp(mol_for_fit{imol},'HDO')
%                     eVMR = eVMR*3.106930e-4;
%                 end
                inp_interp.Pq = Csublayer(isublayer,1);% sublayer pressure, hPa
                inp_interp.Tq = Csublayer(isublayer,3);% sublayer temperature, K
                
                eZ = Csublayer(isublayer,2)*1e5;% sublayer thickness, km->cm
                eN = eVMR*(inp_interp.Pq*100)/kB/inp_interp.Tq*1e-6; % sublayer number density, mol/cm3
                eTau = eN*F_interp_absco(inp_interp)*eZ; % sublayer optical depth
                eTau(isnan(eTau)) = 0;
                count_interp = count_interp+1;
                dTau_layer = dTau_layer+eTau;
            end
            % convolve and interpolate sum of sublayers into each layer, at
            % common grid
            Tau_layer(ilevel,:) = ...
                F_conv_interp(inp_interp.Wq,dTau_layer,dgrd_fwhm,common_grid);
        end
        tmp_struct.(mol_for_fit{imol}).Tau_sum = sum(Tau_layer);
        
        % calculate the absorption cross-section at the sfc layer
        inp_interp.Tq = mean([C(end,3) C_0(end,3)]);
        inp_interp.Pq = inp.surface_layer_P;
        tmp_struct.(mol_for_fit{imol}).surface_layer_sigma = ...
            F_conv_interp(inp_interp.Wq,...
            F_interp_absco(inp_interp),dgrd_fwhm,common_grid);
        
        tmp_struct.(mol_for_fit{imol}).surface_layer_top_P = ...
            inp.lowest_possible_Psurf;
        
    end
    window_list(iwin).tau_struct = tmp_struct;
    
    % add O2 CIA to the first window
    if iwin == 1
        datadir = inp.CIA_dir;
        if isfield(inp,'which_CIA')
            which_CIA = inp.which_CIA;
        else
            which_CIA = 'gfit';
        end
        
        if ~iscell(which_CIA)
            which_CIA = {which_CIA};
        end
        
        % load O2 CIA at 1.27 um band
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
            C_mate = cell(2,1);
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
                C_mate(count) = textscan(fid,'%f%f','Delimiter',' ','MultipleDelimsAsOne',1,...
                    'collectoutput',1);
            end
            fclose(fid);
            ciaT_mate = ciaT;
            vlow_mate = C_mate{1}(:,1);
            C_mate_old = C_mate;
            c_length = length(C_mate_old);
            C_mate = repmat(C_mate_old{1}(:,2),[1,length(C_mate_old)]);
            for imate = 2:c_length
                C_mate(:,imate) = interp1(C_mate_old{imate}(:,1),C_mate_old{imate}(:,2),vlow_mate,'linear','extrap');
            end
            [xgrid_mate,ygrid_mate] = meshgrid(vlow_mate,[253 273 296]);
        end
        if sum(ismember(which_CIA,'gfit'))
            dtau_cntm_gfit = zeros(nlayer,length(common_grid),'single');
        end
        if sum(ismember(which_CIA,'sao'))
            dtau_cntm_sao = zeros(nlayer,length(common_grid),'single');
        end
        if sum(ismember(which_CIA,'mate'))
            dtau_cntm_mate = zeros(nlayer,length(common_grid),'single');
        end
        for ilayer = 1:nlayer
            MR_O2 = C(ilayer,7); % has to be profile of O2, 0.2095
            MR_N2 = 0.78;
            % number density of air in this layer in molec/cm3
            N_ilayer = C(ilayer,1)*100/1e6/kB/C(ilayer,3);
            if sum(ismember(which_CIA,'gfit'))
                [~, fcia_temp] = ...
                    F_CIA2Spec(gfit_fcia,C(ilayer,3),C(ilayer,1),vStart,vEnd,1,v_grid_gfit);
                [~, scia_temp] = ...
                    F_CIA2Spec(gfit_scia,C(ilayer,3),C(ilayer,1),vStart,vEnd,1,v_grid_gfit);
                fcia_low = accumarray(bin(:),fcia_temp,[],@mean);
                scia_low = accumarray(bin(:),scia_temp,[],@mean);
                acia_low = MR_O2*scia_low + MR_N2*fcia_low;
                dtau_cntm_gfit(ilayer,:) = N_ilayer^2*MR_O2*1e5*(C(ilayer,2)-C(ilayer+1,2))*...
                    interp1(0.5*(vlow_gfit(1:end-1)+vlow_gfit(2:end)),acia_low,common_grid);
            end
            if sum(ismember(which_CIA,'sao'))
                if C(ilayer,3) >= max(cia_temp_vec)
                    acia_low = C_sao(:,2);
                elseif C(ilayer,3) <= min(cia_temp_vec)
                    acia_low = C_sao(:,end);
                else
                    acia_low = interp2(xgrid,ygrid,C_sao(:,2:end)',vlow_sao,C(ilayer,3));
                end
                dtau_cntm_sao(ilayer,:) = N_ilayer^2*MR_O2*1e5*(C(ilayer,2)-C(ilayer+1,2))*...
                    interp1(vlow_sao,acia_low,common_grid);
            end
            if sum(ismember(which_CIA,'mate'))
                if C(ilayer,3) >= 296
                    acia_low = C_mate(:,3);
                elseif C(ilayer,3) <= 253
                    acia_low = C_mate(:,1);
                else
                    acia_low = interp2(xgrid_mate,ygrid_mate,C_mate',vlow_mate,C(ilayer,3));
                end
                dtau_cntm_mate(ilayer,:) = N_ilayer^2*MR_O2*1e5*(C(ilayer,2)-C(ilayer+1,2))*...
                    interp1(vlow_mate,acia_low,common_grid);
            end
        end
        
        if sum(ismember(which_CIA,'gfit'))
            window_list(iwin).tau_struct.O2.CIA_GFIT = nansum(dtau_cntm_gfit);
        end
        if sum(ismember(which_CIA,'sao'))
            window_list(iwin).tau_struct.O2.CIA_SAO = nansum(dtau_cntm_sao);
        end
        if sum(ismember(which_CIA,'mate'))
            window_list(iwin).tau_struct.O2.CIA_Mate = nansum(dtau_cntm_mate);
        end
    end
    disp(['Finished ',target_gas,' window ',num2str(windowID),' at ',datestr(now)])
end
%%
% iwin = 2;
% plot(window_list(iwin).common_grid,window_list(iwin).tau_struct.CH4.surface_layer_sigma)