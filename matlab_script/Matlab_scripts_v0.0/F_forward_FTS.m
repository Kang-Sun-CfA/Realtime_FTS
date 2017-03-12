function s_fit = F_forward_FTS(coeff,input)
% forward function for EM27. Written by Kang Sun on 2016/07/17
ils = input.ils;
w1 = input.w1;
ss = input.ss;
s1 = input.s1;
w2 = input.w2;
nvar = 1;
if ismember('cntm_level',input.varlist)
    cntm_level = coeff(nvar);nvar = nvar+1;
else
    cntm_level = 1;
end
if ismember('cntm_tilt',input.varlist)
    cntm_tilt = coeff(nvar);nvar = nvar+1;
else
    cntm_tilt = 0;
end
if ismember('gas_shift',input.varlist)
    gas_shift = coeff(nvar);nvar = nvar+1;
else
    gas_shift = 0;
end
if ismember('solar_shift',input.varlist)
    solar_shift = coeff(nvar);nvar = nvar+1;
else
    solar_shift = 0;
end
if ismember('ZLO',input.varlist)
    ZLO = coeff(nvar);nvar = nvar+1;
else
    ZLO = 0;
end
if ismember('solar_scaling',input.varlist)
    solar_scaling = coeff(nvar);nvar = nvar+1;
else
    solar_scaling = 1;
end
if ismember('CIA_scaling',input.varlist)
    CIA_scaling = coeff(nvar);nvar = nvar+1;
else
    CIA_scaling = 1;
end
ss_over = solar_scaling*conv(ss,ils/sum(ils),'same');
s1_all = zeros(1,size(s1,2));
for imol = 1:size(s1,1)
s1_all = s1_all+s1(imol,:)*coeff(nvar-1+imol);
end
s1_over = conv(s1_all,ils/sum(ils),'same');

s_fit = ZLO+polyval([cntm_tilt,cntm_level],w2-mean(w2))...
    .*interp1(w1,ss_over,w2+solar_shift,'linear','extrap').*...
    exp(-interp1(w1,s1_over,w2+gas_shift,'linear','extrap'));
if isfield(input,'s_CIA')
%     s_CIA = input.s_CIA;
%     s_CIA_over = conv(s_CIA(1,:)*CIA_scaling,ils/sum(ils),'same');
    s_CIA_over = input.s_CIA*CIA_scaling;
    s_fit = s_fit.*...
        exp(-interp1(w1,s_CIA_over,w2+gas_shift,'linear','extrap'));
end
s_fit = s_fit(:);
    

