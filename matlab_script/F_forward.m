function s_fit = F_forward(coeff,input_nlinfit)
% forward function for EM27. Written by Kang Sun on 2016/07/17
% rewritten from F_forward_FTS.m as F_forward by Kang Sun on 2017/01/25

% must-have inputs
ils = input_nlinfit.ils;
w1 = input_nlinfit.w1;
ss = input_nlinfit.ss;
s1 = input_nlinfit.s1;
w2 = input_nlinfit.w2;

% add fitting boundaries
if isfield(input_nlinfit,'coeff_bound')
    tmp = input_nlinfit.coeff_bound;
    for ib = 1:size(tmp,2)
        if coeff(ib) < tmp(1,ib)
            coeff(ib) = tmp(1,ib);
        elseif coeff(ib) > tmp(2,ib);
            coeff(ib) = tmp(2,ib);
        end
    end
end

% optional inputs

count = 1;

if ~isfield(input_nlinfit,'shift_gas')
    shift_gas = coeff(count);
    count = count+1;
else
    shift_gas = input_nlinfit.shift_gas;
end

if ~isfield(input_nlinfit,'shift_sun')
    shift_sun = coeff(count);
    count = count+1;
else
    shift_sun = input_nlinfit.shift_sun;
end

if ~isfield(input_nlinfit,'shift_together')
    shift_together = coeff(count);
    count = count+1;
else
    shift_together = input_nlinfit.shift_together;
end

if ~isfield(input_nlinfit,'scaling')
    scaling = coeff(count);
    count = count+1;
else
    scaling = input_nlinfit.scaling;
end

if ~isfield(input_nlinfit,'tilt')
    tilt = coeff(count);
    count = count+1;
else
    tilt = input_nlinfit.tilt;
end

if ~isfield(input_nlinfit,'zlo')
    zlo = coeff(count);
    count = count+1;
else
    zlo = input_nlinfit.zlo;
end

s1_all = zeros(1,size(s1,2));
for imol = 1:size(s1,1)
    s1_all = s1_all+s1(imol,:)*coeff(count-1+imol);
end

if shift_sun ~= 0
    ss = interp1(w1,ss,w1+shift_sun,'linear','extrap');
end
if shift_gas ~= 0
    s1_all = interp1(w1,s1_all,w1+shift_gas,'linear','extrap');
end
transmission = conv(ss.*exp(-s1_all),ils/sum(ils),'same');
s_fit = zlo+polyval([tilt, scaling],w2-mean(w2))...
    .*interp1(w1,transmission,w2+shift_together,'linear','extrap');
s_fit = s_fit(:);

% % older version, I0 effect could be an issue

% ss_over = coeff(6)*conv(ss,ils/sum(ils),'same');
% s1_all = zeros(1,size(s1,2));
% for imol = 1:size(s1,1)
% s1_all = s1_all+s1(imol,:)*coeff(6+imol);
% end
% s1_over = conv(s1_all,ils/sum(ils),'same');
% 
% s_fit = coeff(5)+polyval(coeff([2 1]),w2-mean(w2))...
%     .*interp1(w1,ss_over,w2+coeff(4),'linear','extrap').*...
%     exp(-interp1(w1,s1_over,w2+coeff(3),'linear','extrap'));
% if isfield(input_nlinfit,'s_CIA')
%     s_CIA = input_nlinfit.s_CIA;
%     s_CIA_over = conv(s_CIA(1,:)*coeff(end),ils/sum(ils),'same');
%     s_fit = s_fit.*...
%         exp(-interp1(w1,s_CIA_over,w2+coeff(3),'linear','extrap'));
% end
% s_fit = s_fit(:);