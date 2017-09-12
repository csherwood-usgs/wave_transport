function ksw = ksw_func( d50, rh, rl, theta )
% Calculate wave roughess
% ksw = ksw_func( d50, rh, rl, theta )
%
% Input:
% d50 - median grain size (m)
% rh - ripple height (m)
% rh - ripple wavelength (m)
% theta - time-averaged absolute Shields stress ()
%
% Returns:
% ksw - roughness (m)
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

% eqn A.5
ksw = max(d50, d50*(mu+6.*(theta-1.))) + 0.4*rh*rh/rl;
return