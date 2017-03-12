clc

% load input parameters
S_input_parameters;
if ~exist(NAMdatapath,'dir')
    mkdir(NAMdatapath);
end
cd(homepath)
addpath(genpath(functionpath))

if_download_NAM = false;
if_extract_NAM = false;
if_plot_NAM = false;
if_line_by_line = true;
if_test_retrieval = false;
if if_download_NAM
    parfor ihr = 1:length(hour_ahead)
        tic
        webfilename = ['nam.t',hour_run,'z.conusnest.hiresf',...
            sprintf('%0*d',2,hour_ahead(ihr)),'.tm00.grib2'];
        url = ['http://www.ftp.ncep.noaa.gov/data/nccf/com/nam/prod/nam.'...
            year_run,month_run,day_run,'/',webfilename];
        
        outfilename = websave([NAMdatapath,webfilename],url);
        disp(['The file ',webfilename,' is saved at ',NAMdatapath,'.']);
        toc
    end
end
% define map projection and locate the measurement site
mstruct = defaultm('lambertstd');
mstruct.origin = [25 265 0];
mstruct.mapparallels = 25;
mstruct.nparallels = 1;
mstruct.scalefactor = 6371229/1e3;
mstruct.falseeasting = 0;
mstruct.falsenorthing = 0;
mstruct = defaultm(mstruct);

[target_x, target_y] = mfwdtran(mstruct,target_lat,target_lon);

NSA = shaperead('landareas.shp','UseGeoCoords',true,...
    'Selector',{@(name)strcmpi(name,'North and South America'),'Name'});
[NSA_x,NSA_y] = mfwdtran(mstruct,NSA.Lat,NSA.Lon);

if if_extract_NAM
    namP_mat = nan(42,length(hour_ahead));
    namT_mat = namP_mat;
    namH_mat = namP_mat;
    namH2O_mat = namP_mat;
    namP_surf_vec = nan(1,length(hour_ahead));
    namT_surf_vec = namP_surf_vec;
    parfor ihr = 1:length(hour_ahead)
        webfilename = ['nam.t',hour_run,'z.conusnest.hiresf',...
            sprintf('%0*d',2,hour_ahead(ihr)),'.tm00.grib2'];
        try
            namData = F_read_nam([NAMdatapath,webfilename],target_x,target_y);
        catch
            setup_nctoolbox
            namData = F_read_nam([NAMdatapath,webfilename],target_x,target_y);
        end
        namP_mat(:,ihr) = namData.isobaric;
        namT_mat(:,ihr) = namData.namT;
        namH_mat(:,ihr) = namData.namH;
        namH2O_mat(:,ihr) = namData.namH2O;
        
        namP_surf_vec(ihr) = namData.namP_surf;
        namT_surf_vec(ihr) = namData.namT_surf;
        if if_plot_NAM
            close all
            figure('color','w','unit','inch','position',[0 0 7.5 6])
            h = pcolor(namData.x,namData.y,double(namData.Temperature_surface_all));
            set(h,'edgecolor','none')
            hold on; plot(NSA_x,NSA_y,'w','linewidth',1.5)
            plot(target_x,target_y,'-pr','markersize',12,'linewidth',1);
            
            
            hold off;text(target_x-700,target_y+200,'Boston, MA')
            hc = colorbar;set(get(hc,'ylabel'),'string','Surface T [K]')
            caxis([270 315])
            xlabel('Distance [km]')
            ylabel('Distance [km]')
            title(['NAM forecast made at ',month_run,'/',day_run,', ',hour_run,...
                ':00 UTC, for ',sprintf('%0*d',2,hour_ahead(ihr)),' hours ahead']);
            export_fig(['/data/wdocs/kangsun/www-docs/files/temp',num2str(ihr),'.png'],'-r150');
        end
    end
end

if if_plot_NAM
    close all
    plot(namT_mat,namP_mat)
    set(gca,'yscale','log','ydir','rev')
    export_fig(['/data/wdocs/kangsun/www-docs/files/ensemble_profT.png'],'-r150')
    close all
    plot(namH_mat,namP_mat)
    set(gca,'yscale','log','ydir','rev')
    export_fig(['/data/wdocs/kangsun/www-docs/files/ensemble_profH.png'],'-r150')
    close all
    plot(namH2O_mat,namP_mat)
    set(gca,'yscale','log','ydir','rev')
    export_fig(['/data/wdocs/kangsun/www-docs/files/ensemble_profH2O.png'],'-r150')
    save([NAMdatapath,'fcst_prof',year_run,month_run,day_run,hour_run,'.mat'],'hour_ahead','namP_mat','namT_mat','namH_mat','namH2O_mat','namP_surf_vec','namT_surf_vec')
    
    run_UTC = datenum([year_run,'-',month_run,'-',day_run,' ',hour_run,':00:00']);
    gifname = '/data/wdocs/kangsun/www-docs/files/surfT_animation.gif';
    delaytime = .9;
    for t= 1:length(hour_ahead)
        Frame=imread(['/data/wdocs/kangsun/www-docs/files/temp',num2str(t),'.png']);
        forecast_UTC = run_UTC+hour_ahead(t)/24;
        % im = AddTextToImage(Frame,datestr(forecast_UTC),...
        %    [260 900],[0 0 0],'Tahoma',45);
        im = Frame;
        image(im)
        axis off
        [A,map] = rgb2ind(im,256);
        if t == 1;
            imwrite(A,map,gifname,'gif','LoopCount',Inf,'DelayTime',delaytime);
        else
            imwrite(A,map,gifname,'gif','WriteMode','append','DelayTime',delaytime);
        end
    end
end
if if_line_by_line
    load([NAMdatapath,'fcst_prof',year_run,month_run,day_run,hour_run,'.mat'],'hour_ahead','namP_mat','namT_mat','namH_mat','namH2O_mat','namP_surf_vec','namT_surf_vec')
    %% Hack GFIT vmr profiles
    fid = fopen('/data/tempo1/Shared/kangsun/ggg_location/vmrs/gnd/gnd_summer.vmr');
    C = textscan(fid,['%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f',...
        '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f'],...
        'headerlines',5);
    C = cell2mat(C);
    fclose(fid);
    profileinput.C = C;
    %% load solar spectrum
    load([datadir, 'FTS_solar.mat'],'llwaven','lltranswaven')
    
    for ihr = 1:1%length(hour_ahead)
        % assume geometric height is geopotential height, Is that OK?
        % if not OK, need a function to tweak namData.namH
        savepath = ['/data/tempo1/Shared/kangsun/FTS_data/Tau_data/',year_run,...
            month_run,day_run,'/UTC',hour_run,'f',num2str(hour_ahead(ihr)),'/'];
        if ~exist(savepath,'dir')
            mkdir(savepath);
        end
        Z_level = namH_mat(:,ihr);
        dZ = diff(Z_level);
        P_level = namP_mat(:,ihr)/100; % pressure in hPa
        P_layer = P_level(1:end-1)+diff(P_level)/2;
        nlayer = length(P_layer);
        T_level = namT_mat(:,ihr);
        T_layer = interp1(P_level,T_level,P_layer);
        Z_layer = interp1(P_level,Z_level,P_layer);
        profileinput.Z_layer = Z_layer;
        profileinput.cfsH2O = namH2O_mat(:,ihr);
        
        clc
        for iwin = 1:length(window_list)
            target_gas = window_list(iwin).target_gas;
            windowID = window_list(iwin).windowID;
            [vStart, vEnd, wStart, wEnd, molnameinput,mol_for_spec] = ...
                F_define_windows(target_gas,windowID);
            int = llwaven >= vStart & llwaven <= vEnd;
            llwave = llwaven(int);lltrans = lltranswaven(int);
            % use llwave as common_grid
            common_grid = llwave;
            save([savepath,'common_grid_',num2str(vStart),...
                '_',num2str(vEnd),'.mat'],'common_grid')
            for imol = 1:length(mol_for_spec)
                taumol = mol_for_spec{imol};
                output_LBL = ...
                    F_line_by_line_FTS(taumol,...
                    F_ap_profile(taumol,profileinput),P_level,P_layer,T_layer,...
                    Z_level,Z_layer,vStart,vEnd,molparampath,datadir,common_grid);
                save([savepath,'tau_data_',taumol,'_',num2str(vStart),...
                    '_',num2str(vEnd),'.mat'],'output_LBL')
                disp(['Finished forecast hour ',num2str(hour_ahead(ihr)),...
                    ', ',target_gas,' window #',num2str(windowID),...
                    ',',char(10),'molecule ',taumol,' at ',datestr(now),char(10)])
            end
            % save([savepath,'ss_',num2str(vStart),'_',num2str(vEnd),'.mat'],...
            %     'llwave','lltrans')
        end
    end
end
%%
if if_test_retrieval
    load([NAMdatapath,'fcst_prof',year_run,month_run,day_run,hour_run,'.mat'],'hour_ahead','namP_mat','namT_mat','namH_mat','namH2O_mat','namP_surf_vec','namT_surf_vec')
    
    clc
    close all
    ifile = 100;
    filename = [spectrapath,'ha20150807.s0e00a.',sprintf('%04d',ifile)];
    [obs_y,obs_x,params] = ImportOpus(filename,'SampleSpectrum');
    load([datadir, 'FTS_solar.mat'],'llwaven','lltranswaven')
    for ihr = 1:1
        Z_level = namH_mat(:,ihr);
        dZ = diff(Z_level);
        P_level = namP_mat(:,ihr)/100; % pressure in hPa
        P_layer = P_level(1:end-1)+diff(P_level)/2;
        nlayer = length(P_layer);
        T_level = namT_mat(:,ihr);
        T_layer = interp1(P_level,T_level,P_layer);
        Z_layer = interp1(P_level,Z_level,P_layer);
        % assume geometric height is geopotential height, Is that OK?
        % if not OK, need a function to tweak namData.namH
        savepath = ['/data/tempo1/Shared/kangsun/FTS_data/Tau_data/',year_run,...
            month_run,day_run,'/UTC',hour_run,'f',num2str(hour_ahead(ihr)),'/'];
        for iwin = 1:length(window_list)
            target_gas = window_list(iwin).target_gas;
            windowID = window_list(iwin).windowID;
            [vStart, vEnd, wStart, wEnd, mol_for_fit,mol_for_spec] = ...
                F_define_windows(target_gas,windowID);
            load([savepath,'common_grid_',num2str(vStart),...
                '_',num2str(vEnd),'.mat'],'common_grid')
            dtau_mol = nan(length(mol_for_fit),length(common_grid),nlayer);
            % clear the fitting input from the previous window
            clear input
            for imol = 1:length(mol_for_fit)
                taumol = mol_for_fit{imol};
                if strcmpi(taumol,'HDO')
                    load([savepath,'tau_data_','H2O','_',num2str(vStart),...
                        '_',num2str(vEnd),'.mat'],'output_LBL')
                    dtau_mol(imol,:,:) = output_LBL.dtau_layer(2,:,:);
                elseif strcmpi(taumol,'O2')
                    load([savepath,'tau_data_',taumol,'_',num2str(vStart),...
                        '_',num2str(vEnd),'.mat'],'output_LBL')
                    dtau_mol(imol,:,:) = output_LBL.dtau_layer(1,:,:);
                    % static input to forward function
                    input.s_CIA = double(sum(output_LBL.tau_cntm,1));
                else
                    load([savepath,'tau_data_',taumol,'_',num2str(vStart),...
                        '_',num2str(vEnd),'.mat'],'output_LBL')
                    dtau_mol(imol,:,:) = output_LBL.dtau_layer(1,:,:);
                end
            end
            int = llwaven >= vStart & llwaven <= vEnd;
            llwave = llwaven(int);lltrans = lltranswaven(int);
            int = obs_x >= wStart & obs_x <= wEnd;
            w2 = double(obs_x(int));
            s2 = double(obs_y(int));
            % use ideal, boxcar ILS.
            vCenter = mean([vStart vEnd]);
            ilsextent = 50; % wavenumber
            ilsx = -ilsextent:median(diff(common_grid)):ilsextent;
            OPD = 1.8; % maximum opd in cm
            FOV = 2.96e-3; % field of view in rad
            ils_sync = 2*OPD*sinc(2*ilsx*OPD);
            ils_rect = zeros(size(ilsx));
            ils_rect(ilsx >= -.5*vCenter*FOV^2 & ilsx <= 0) = 2/(vCenter*FOV^2);
            modulation_loss = 0.025;
            ils_align = sinc(ilsx*OPD/modulation_loss).^2;
            ilsy = conv(ils_sync,ils_rect,'same');
            % static input to forward function
            input.ils = ilsy;
            input.w2 = w2;
            input.w1 = llwave;
            input.ss = interp1(llwave,lltrans,common_grid);
            input.s1 = squeeze(sum(dtau_mol,3));

            
            % state vector:
            % 1 continuum level
            % 2 continuum tilt
            % 3 gas wave number shift
            % 4 solar wavenumber shift
            % 5 ZLO
            % 6 solar scaling
            coeff0 = [prctile(s2,90) 0 0 0 0 1];
            % 7+ add target gas scaling
            coeff0 = [coeff0 ones(1,length(mol_for_fit))];
            % add O2 CIA scaling, if there is O2
            if ismember('O2',mol_for_fit)
                coeff0 = [coeff0 1];
            end
            tic
            [coeff, R] = nlinfit(input,s2,@F_forward_FTS,coeff0);
            fittime = toc;
            disp(['Non linear fitting took ',num2str(fittime),' s'])
            
            w1 = common_grid;
            figure('unit','inch','position',[-10 1 6.5 9],'color','w')
            ax = cell(3,1);
            ax{1} = subplot(3,1,1);
            % h = plot(w2,s2,'k',w2,F_forward_FTS(coeff0,input),w2,F_forward_FTS(coeff,input));
            
            h1 = plot(w2,s2,'k',w2,F_forward_FTS(coeff,input),'r','linewidth',1);
            hleg = legend(h1,'Observation','Fitting');
            set(hleg,'location','southwest','box','off')
            ax{2} = subplot(3,1,2);
            if strcmp(target_gas,'O2')
                h2 = plot(w1,1-exp(-input.s1),w1,input.s_CIA(1,:),w1,input.ss,'linewidth',1);
                hleg = legend(h2,{mol_for_fit{:},'O2-O2 CIA','Solar transmission'});
            else
                h2 = plot(w1,1-exp(-input.s1),w1,input.ss,'linewidth',1);
                hleg = legend(h2,{mol_for_fit{:},'Solar transmission'});
            end
            set(hleg,'location','west','box','off')
            ax{3} = subplot(3,1,3);
            plot(w2,R/coeff(1)*100,'k','linewidth',1)
            hleg = legend('Fitting residual [% of continnum]');
            set(hleg,'location','southwest','box','off')
            for iplot = 1:3
                set(ax{iplot},'xlim',[wStart wEnd],'linewidth',1)
                if iplot < 3
                    set(ax{iplot},'xcolor','none')
                else
                    xlabel('Wave number [cm^{-1}]')
                end
            end
        end
    end
end