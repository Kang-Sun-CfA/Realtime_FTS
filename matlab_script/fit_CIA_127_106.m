% overhaul the fit CIA code to use absco tables
% updated by Kang Sun from fit_CIA_106.m and F_absco2tau.m on 2017/10/07

clear;clc;close all
% system separator, "/" for mac and linux, "\" for windows
sfs = filesep;
inp = [];

%!!!!!!!!!!!!!!!! local dir saving ABSCO tables !!!!!!!!!!!!!!!!!!!!!!
inp.absco_dir = '/data/tempo1/Shared/kangsun/FTS_data/HITRAN/';
%!!!!!!!!!!!!!!!! YOU HAVE TO SPECIFY THIS !!!!!!!!!!!!!!!!!!!!!!!!!!!

inp.CIA_dir = ['..',sfs,'spectroscopy',sfs];
inp.which_CIA = {'sao','gfit','mate'};
inp.profile_fn = 'D:\Research_CfA\FTS\20160624\hb20160624.map';
inp.lowest_possible_Psurf = 980;
inp.nsublayer = 1;
inp.common_grid_resolution = 0.05;
inp.surface_layer_P = 990;

% end of input

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

%% load profile
fid = fopen(inp.profile_fn);
C_0 = cell2mat(textscan(fid,repmat('%f',[1,12]),'headerlines',11,'delimiter',',','multipledelimsasone',1));
% I insist that the first column should be pressure and the second column
% be altitude. The third one better to be temperature. 
C0 = C_0;
C0(:,[1 2 3 4 5 6]) = C_0(:,[3 1 2 5 7 10]);
% I don't like definition of co2 and ch4. Change them to vmr
C0(:,5) = C0(:,5)*1e-6;
C0(:,6) = C0(:,6)*1e-9;
% add O2 mixing ratio
C0(:,7) = ones(size(C0,1),1)*0.2095;
C0 = C0(:,1:7);
C_0 = C0;

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
%% make some plots
if ispc
load('C:\Users\Kang Sun\Documents\GitHub\Realtime_FTS\spectroscopy\ptgrid.mat')
load('C:\Users\Kang Sun\Documents\GitHub\Realtime_FTS\spectroscopy\abscoptgrid.mat')

fid = fopen('D:\Research_CfA\FTS\20170217\hb20170217.map');
C_0 = cell2mat(textscan(fid,repmat('%f',[1,12]),'headerlines',11,'delimiter',',','multipledelimsasone',1));
% I insist that the first column should be pressure and the second column
% be altitude. The third one better to be temperature. 
C0 = C_0;
C0(:,[1 2 3 4 5 6]) = C_0(:,[3 1 2 5 7 10]);
% I don't like definition of co2 and ch4. Change them to vmr
C0(:,5) = C0(:,5)*1e-6;
C0(:,6) = C0(:,6)*1e-9;
% add O2 mixing ratio
C0(:,7) = ones(size(C0,1),1)*0.2095;
C0 = C0(:,1:7);
C_0 = C0;

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

C_0217 = cat(1,C_0(C_0(:,1)-inp.lowest_possible_Psurf < 0,:),C_low);

close all
set(0,'defaultaxesfontsize',12)
figure('unit','inch','color','w','position',[15 1 14 6])
subplot(1,2,1)
semilogy(tempgrid_absco,presgrid_absco,'.k','markersize',12);axis ij
hold on
hw = plot(C_0217(:,3),C_0217(:,1)*100,'b','linewidth',3);
hs = plot(C(:,3),C(:,1)*100,'r','linewidth',3);
hold off
set(gca,'xlim',[150 350],'ylim',[1 1.2e5],'linewidth',1,'box','off')
ylabel('Pressure [Pa]')
xlabel('Temperature [K]')
hleg = legend([hw hs],'Winter (Feb 17, 2017)','Summer (Jun 24, 2016)');
set(hleg,'box','off')
title('OCO-2 ABSCO P/T grid','fontsize',16)

subplot(1,2,2)
semilogy(tempgrid,presgrid,'.k','markersize',12);axis ij
hold on
hw = plot(C_0217(:,3),C_0217(:,1)*100,'b','linewidth',3);
hs = plot(C(:,3),C(:,1)*100,'r','linewidth',3);
hold off
set(gca,'xlim',[150 350],'ylim',[1 1.2e5],'linewidth',1,'box','off')
xlabel('Temperature [K]')
title('Modified P/T grid','fontsize',16)
%%
close all
set(0,'defaultaxesfontsize',12)
ivar = 3;
switch ivar
    case 2
Xlim = [-0.3 2.1];
offset = 0.1;
    case 3
        Xlim = [270 310];
        offset = 1.5;
end
Ylim = [780 1050];
figure('unit','inch','color','w','position',[-15 1 12 6])
subplot(1,3,1)
hold on
plot(C_0(:,ivar),C_0(:,1),'-O','linewidth',2);axis ij
count = 0;
for i = size(C_0,1):-1:size(C_0,1)-3
    plot([Xlim(1) C_0(i,ivar)],[C_0(i,1),C_0(i,1)],':k')
    text(C_0(i,ivar)+offset,C_0(i,1),['Level ',num2str(count)],...
        'horizontalalignment','left','fontsize',16);
    count = count+1;
end
plot(C_low(ivar),C_low(1),'rp','linewidth',2,'markersize',16)
text(C_low(ivar)+offset,C_low(1)+0.,['Lowest possible',char(10), 'surface pressure'],...
        'horizontalalignment','left','fontsize',14);

hold off
ylim(Ylim)
xlim(Xlim)
title('A priori profiles','fontsize',16)
set(gca,'linewidth',1)
ylabel('Pressure [hPa]')
subplot(1,3,2)
hold on
plot(C(:,ivar),C(:,1),'-O','linewidth',2);axis ij
count = 0;
for i = size(C,1):-1:size(C,1)-3
    plot([Xlim(1) C(i,ivar)],[C(i,1),C(i,1)],':k')
    text(C(i,ivar)+offset,C(i,1),['Level ',num2str(count)],...
        'horizontalalignment','left','fontsize',16);
    count = count+1;
end
plot(C_low(ivar),C_low(1),'rp','linewidth',2,'markersize',16)
text(C_low(ivar)-offset,C_low(1)+0.,['Lowest possible',char(10), 'surface pressure'],...
        'horizontalalignment','right','fontsize',14);

hold off
ylim(Ylim)
xlim(Xlim)
title('Truncated profiles','fontsize',16)
set(gca,'linewidth',1,'ycolor','w','ytick',[])
xlabel('Temperature [km], but can be VMRs')
% ylabel('Relative altitude [km]')

subplot(1,3,3)
hold on
hold on
plot(C(:,ivar),C(:,1),'-O','linewidth',2);axis ij
count = 0;
for i = size(C,1):-1:size(C,1)-3
    plot([Xlim(1) C(i,ivar)],[C(i,1),C(i,1)],':k')
    text(C(i,ivar)+offset,C(i,1),['Level ',num2str(count)],...
        'horizontalalignment','left','fontsize',16);
    count = count+1;
end

P_surface = 999;
plot(C(end,ivar),P_surface,'k^','linewidth',2,'markersize',10);

text(C(end,ivar)+offset,P_surface+10,['Real-time',char(10),'surface pressure'],...
        'horizontalalignment','left','fontsize',14);
    plot([Xlim(1) C(end,ivar)],[P_surface,P_surface],':k')
ha = area([Xlim(1) C(end,ivar)],[C(end,1),C(end,1);P_surface-C(end,1),P_surface-C(end,1)]');
set(ha(1),'facecolor','none','edgecolor','none')
set(ha(2),'facecolor','c','edgecolor','none')
uistack(ha,'bottom')
hold off
ylim(Ylim)
xlim(Xlim)
title('Truncated profiles+adjustable surface layer','fontsize',16)
set(gca,'linewidth',1,'ycolor','w','ytick',[])
end
%% %%%% Define retrieval windows
field1 = 'target_gas';
field2 = 'windowID';
field3 = 'common_grid';
field4 = 'tau_struct';
field5 = 'vRange';
field6 = 'wRange';
value1 = {'O2','O2'};
value2 = {0,1};
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
        % sfc layer pressure in hPa
        tmp_struct.(mol_for_fit{imol}).surface_layer_top_P = ...
            inp.lowest_possible_Psurf;
        % sfc layer VMR
        tmp_struct.(mol_for_fit{imol}).surface_layer_VMR = C(end,prof_index);
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