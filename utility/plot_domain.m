clear;clc
cd /data/tempo1/Shared/kangsun/run_WRF/WPS
addpath('~/matlab functions/export_fig')
addpath('~/matlab functions/')
plotdir = '/data/tempo1/Shared/kangsun/run_WRF/figures/buffalo/';
%% load state and lake shape files
Llake     = shaperead('/data/tempo1/Shared/kangsun/run_WRF/shapefiles/ne_10m_lakes/ne_10m_lakes.shp');
lakelist = [];
left = -90;right = -74;down = 39.5;up = 46;
for i = 1:length(Llake)
    tmp = Llake(i).BoundingBox;
    if ((tmp(1,1) > left && tmp(1,1) < right) || ...
            (tmp(2,1) > left && tmp(2,1) < right)) && ...
            ((tmp(2,1) > down && tmp(2,1) < up) || ...
            (tmp(2,2) > down && tmp(2,2) < up))
        lakelist = [lakelist i];
    end
end

statelist = [8, 55, 32, 24 46];
S         = shaperead('/data/tempo1/Shared/kangsun/run_WRF/shapefiles/cb_2015_us_state_500k/cb_2015_us_state_500k.shp');
%% load lat, lon, and some 2-D field
domainN = 1:3;
wrfstruct = [];
for i = domainN
    fn = ['geo_em.d0',num2str(domainN(i)),'.nc'];
    wrfout = F_ncread_all(fn);
    lonname = ['wlon_',num2str(domainN(i))];
    latname = ['wlat_',num2str(domainN(i))];
    zname = ['z_',num2str(domainN(i))];
    wrfstruct.(lonname) = wrfout.XLONG_M.data;
    wrfstruct.(latname) = wrfout.XLAT_M.data;
    wrfstruct.(zname) = wrfout.HGT_M.data;
end
%% plot and save
close all
figure('color','w','unit','inch','position',[1 0 11 7.5])
hold on
CC = lines(6);
CC = CC([1 2 4],:);
% lon, lat, and color limit
Xlim = [-93 -65];Ylim = [35.2 52];Clim = [0 1000];
% plot 2-D field, terrain height in this case
for i = domainN
    lonname = ['wlon_',num2str(domainN(i))];
    latname = ['wlat_',num2str(domainN(i))];
    zname = ['z_',num2str(domainN(i))];
    wlon = wrfstruct.(lonname);
    wlat = wrfstruct.(latname);
    wz = wrfstruct.(zname);
    hp = pcolor(wlon,wlat,double(wz));
    set(hp,'edgecolor','none')
    caxis(Clim)
end
hc = colorbar('south');
% set(hc,'ytick',get(hc,'ylim'),'yticklabel',{'0','max'})
set(hc,'position',[0.35 0.07 0.35 0.02])
set(get(hc,'ylabel'),'string','Terrain height [m]')
% plot state lines and lakes
hold on
for istate = 1:length(S)
    plot(S(istate).X,S(istate).Y,'color','k')
end

for ilake = lakelist
    plot(Llake(ilake).X,Llake(ilake).Y,'color','k')
end
% plot WRF domain boundaries
for i = domainN
    lonname = ['wlon_',num2str(domainN(i))];
    latname = ['wlat_',num2str(domainN(i))];
    zname = ['z_',num2str(domainN(i))];
    wlon = wrfstruct.(lonname);
    wlat = wrfstruct.(latname);
    wz = wrfstruct.(zname);
    h = plot(wlon(1,:),wlat(1,:),'.',wlon(end,:),wlat(end,:),'.',...
        wlon(:,1),wlat(:,1),'.',wlon(:,end),wlat(:,end),'.');
    set(h,'marker','none','linestyle','-','linewidth',2.5,'color','w')
    
end
set(gca,'xlim',Xlim,'ylim',Ylim)
axis off
export_fig([plotdir,'buffalo_nest_domain.png'],'-r150')
%% the following code is random trial
close all
i = 13;
zdata = double(squeeze(wrfout.URB_PARAM.data(:,:,i)));
zdata = double(squeeze(wrfout.LANDUSEF.data(:,:,i)));
zdata = double(wrfout.LU_INDEX.data);
figure('color','w','unit','inch','position',[-10 1 10 7],'visible','on')

h = pcolor(wlon,wlat,zdata);colormap(lines)
set(h,'edgecolor','none')
shading interp
hc = colorbar;
hold on
for istate = statelist
    plot(S(istate).X,S(istate).Y,'color',CC(3,:))
end

for ilake = lakelist
    plot(Llake(ilake).X,Llake(ilake).Y,'color',CC(1,:))
end
%%
tmp = size(wlon);
sn = tmp(2);we=tmp(1);
close
plot(wlon(1,:),wlat(1,:),'.',wlon(end,:),wlat(end,:),'.',...
    wlon(:,1),wlat(:,1),'.',wlon(:,end),wlat(:,end),'.')
hold on
i = 19;j= 57;
plot(wlon(i,j),wlat(i,j),'*')
