% cd c:\users\Kang' Sun'\Dropbox\Code\SAO\
% filename = 'D:\Research_CfA\FTS\hb_20160624\hb20160624s0e00a.0310';
cd C:\users\ksun\Dropbox\code\SAO
clc
config = [];
config.apokind = 1;
% config.refraction = 1;
% config.pathpt = 'C:\\FTS\\em27spectra\\160624\\pt';
config.localt = 0.22;
config.loclat = 39.766667;
config.loclon = -86.15;
ifgdir = 'C:\FTS\em27spectra\160624\';
%%
clc
fileN = 1200;%856;
for ifile = 1:length(fileN)
    filename = [ifgdir,'hb20160624s0e00a.',num2str(fileN(ifile))];
    tic
    calmat = F_calmat(filename,config);
    toc
end
%%
clc
input = [];
input.filename = filename;
input.sim_wn = 1.579792480468750e+04/2;
input.phase_error = 0*5.412e-1;
input.apokind = 7;
input.zerofillingfactor = 5;
output = F_calmat_ils(input);
%
nmax = output.nmax;
%%
h = plot(output.ilsx,output.specre/nmax,output.ilsx,output.specim/nmax,...
    output.ilsx,output.ilsy/nmax,'.k');
set(h([1,2]),'linewidth',1)
set(h(3),'linestyle','-')
legend(h,'specre','specim','specrepc')
xlim([input.sim_wn-5,input.sim_wn+5])
%%
set(gca,'yscale','log')
%%
int = calmat.x > 5000 & calmat.x < 13000;
spec7 = calmat.spec(int);
xx = calmat.x(int);
close all;clc
figure('unit','inch','position',[0 1 12 6])
subplot(2,2,1)
plot(1e7./xx,spec7)
xlim([820 934])
subplot(2,2,2)
plot(1e7./xx,spec7)
xlim([980 1110])
subplot(2,2,3:4)
plot(1e7./xx,spec7)
xlim([1155 1330])
xlabel('Wavelength [nm]')
tstr = ['EM27 at UTC ',datestr(calmat.tutc),'; air mass = ',...
    num2str(calmat.am,'%0.3f'),...
    '; water vapor column = 25.88 mm'];
ha = annotation('textbox',[0.3 0.93 0.6 0.06],'string',tstr);
set(ha,'fontsize',12,'edgecolor','none')
%%
addpath('c:\Users\ksun\Dropbox\matlab functions\export_fig\')
export_fig('FTS_spec_1micron.pdf')
%%
tic
calmat = F_calmat(filename,config);
toc
spec1 = calmat.spec(calmat.x > 7700 & calmat.x < 8100);

config.apokind = 7;
% config.refraction = 0;
tic
calmat = F_calmat(filename,config);
toc
spec7 = calmat.spec(calmat.x > 7700 & calmat.x < 8100);
xx = calmat.x(calmat.x > 7700 & calmat.x < 8100);
plot(xx,spec1,xx,spec7)
%%
config.apokind = 7;
% config.refraction = 0;
config.upperlim = 8100;config.lowerlim = 7700;
tic
calmat = F_calmat(filename,config);
toc
int = calmat.x > 5000 & calmat.x < 12000;
xx = 1e7./calmat.x(int);yy = calmat.spec(int);
plot(xx,yy)
