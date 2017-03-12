function s1_low = F_conv_interp(w1,s1,fwhm,common_grid)
% This function convolves s1 with a Gaussian fwhm, resample it to
% common_grid

% Made by Kang Sun on 2016/08/02

slit = fwhm/1.66511;% half width at 1e

dw0 = median(diff(w1));
ndx = ceil(slit*2.7/dw0);
xx = (0:ndx*2)*dw0-ndx*dw0;
kernel = exp(-(xx/slit).^2);
kernel = kernel/sum(kernel);
s1_over = conv(s1, kernel, 'same');
s1_low = interp1(w1,s1_over,common_grid);