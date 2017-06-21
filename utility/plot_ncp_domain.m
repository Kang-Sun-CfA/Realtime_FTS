cd /data/tempo1/Shared/kangsun/WRF/20160930_nc
addpath('~/matlab functions/export_fig')
addpath('~/matlab functions/')
plotdir = '/data/tempo1/Shared/kangsun/run_WRF/figures/buffalo/';
wrfdomain = 1:3;
CC = lines(6);
S = shaperead('/home/kangsun/methane_NCP/China provinces/China provinces.shp');

%% load wrf data
timestart = datenum('Jun-19-2013 8:00:00');
timeend   = datenum('Jun-19-2013 8:00:00');
timestep  = 1;
nfile     = round((timeend-timestart)*24/timestep)+1;
wrfdomain = 1:3;
wrfstruct = [];
close all
figure('color','w','unit','inch','position',[1 0 11 6.5])
hold on
Clim = [0 2200];
for idomain = wrfdomain
ifile = 1;tmp = datevec(timestart+(ifile-1)*timestep/24);
fn = ['wrfout_d0',num2str(idomain),'_',num2str(tmp(1)),'-',...
    num2str(tmp(2),'%02d'),'-',num2str(tmp(3),'%02d'),'_',...
    num2str(tmp(4),'%02d'),':',num2str(tmp(5),'%02d'),':',num2str(tmp(6),'%02d')];
wrfout = F_ncread_all(fn,{'XLONG','XLAT','LU_INDEX','HGT'});
lonname = ['wlon_',num2str(idomain)];
    latname = ['wlat_',num2str(idomain)];
    zname = ['z_',num2str(idomain)];
wlon = wrfout.XLONG.data;
wlat = wrfout.XLAT.data;
wLU = wrfout.LU_INDEX.data;
wHGT = wrfout.HGT.data;
wrfstruct.(lonname) = wlon;
wrfstruct.(latname) = wlat;
hp = pcolor(wlon,wlat,double(wHGT));
    set(hp,'edgecolor','none')
    caxis(Clim)
end
for istate = 1:length(S)
    plot(S(istate).X,S(istate).Y,'color','k')
end
for idomain = wrfdomain
lonname = ['wlon_',num2str(idomain)];
    latname = ['wlat_',num2str(idomain)];
    wlon = wrfstruct.(lonname);
    wlat = wrfstruct.(latname);

    h = plot(wlon(1,:),wlat(1,:),'.',wlon(end,:),wlat(end,:),'.',...
        wlon(:,1),wlat(:,1),'.',wlon(:,end),wlat(:,end),'.');
    set(h,'marker','none','linestyle','-','linewidth',2.5,'color','w')
end
set(gca,'xlim',[109 126],'ylim',[34.3 42])
axis off
hc = colorbar('south');
% set(hc,'ytick',get(hc,'ylim'),'yticklabel',{'0','max'})
set(hc,'position',[0.35 0.07 0.35 0.02])
set(get(hc,'ylabel'),'string','Terrain height [m]')
export_fig([plotdir,'ncp_nest_domain.png'],'-r150')

