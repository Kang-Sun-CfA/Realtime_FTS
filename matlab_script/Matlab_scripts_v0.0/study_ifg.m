cd D:\Research_CfA\FTS\Matlab_scripts
filename = 'D:\Research_CfA\FTS\hb_20160624\hb20160624s0e00a.0001';
[spectrum,v,vparams] = ImportOpus(filename,'SampleSpectrum');
[ifg,t,tparams] = ImportOpus(filename,'SampleInterferogram');
%%
plot(ifg(1:100:end))
