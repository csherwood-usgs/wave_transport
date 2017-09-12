function ksd = ksd_func( d50, d90, rh, rl, theta )
% Calculate current-related bed roughess
% ksd = ksd_func( d50, d90, rh, rl, theta )
% zo = ksd/30
%
% Input:
% d50 - median grain size (m)
% d90 - 90th percentile grain size (m)
% rh - ripple height (m)
% rh - ripple wavelength (m)
% theta - time-averaged absolute Shields stress
%
% Returns:
% ksd - roughness (m)
% 
% Based on Ribberink (1998) in van der A et al. (2013), Appendix A.
rh = max(rh,d50);
rl = max(rl,d50); % avoid divide-by-zero
% Eqn. A.2
mu = 6.;
d50mm = d50*1.e3;
if( d50mm > 0.15 )
   mu = 6. - 5.*(d50mm-0.15)/(0.2-0.15);
   if( d50mm >= 0.2 )
      mu = 1.;
   end
end
% eqn A.1
ksd = max(3.*d90, d50*(mu+6.*(theta-1.))) + 0.4*rh*rh/rl;
return