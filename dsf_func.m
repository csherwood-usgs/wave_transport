function dsf = dsf_func(d50, theta )
% Calculate sheet flow thickess
% dsf = dsf_func(d50, theta )
%
% Input:
% d50 - median grain size (m)
% theta - maximum (crest or trough) Shields paramter
% Returns:
% dsf - thickness of sheet flow layer (m)
% 
% Based on Janssen (1999) in van der A et al. (2013), Appendix C.
% See also Dohmen-Janssen et al. (2001), JGR 106(C11):27,103-27,115, Eqn. 6 & 7
d50mm = 1.e3*d50;
deln = 25.*theta;
if (d50mm > 0.15 )
   deln = deln-12.*(d50mm-0.15)/(0.2-0.15);
   if (d50mm >= 0.2 )
      deln = 13.*theta;
   end
end
dsf=max( deln*d50, d50 ); % unstated, but assume dsf = d50 when theta = 0
return