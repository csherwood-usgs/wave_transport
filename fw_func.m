function fw = fw_func(ahat,ksw)
% Calculate wave friction factor
% fw = fw_func(ahat,ksw)
fw = 0.3;
aksw = ahat/ksw;
if aksw > 1.587
   fw = 0.00251*exp(5.21*aksw.^(-0.19)); % Eqn. A4, A5
end
return