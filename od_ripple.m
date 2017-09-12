function [Hoa,Loa] = od_ripple( d50, Psi )
% Calculate 0'Donoghue et al. (2006) ripple geometry
% odr = od_ripple( d50, Psi )
%
% Input:
% d50 - median grain size (m)
% Psi - mobility number, max of crest or trough flow
% Returns:
% Hoa - ripple height normalized by orbital amplitude A
% Hoa - ripple length normalized by orbital amplitude A
%
% van der A et al. (2013) Appendix B.
d50mm = 1.e3*d50; %convert from m to mm
mL = 0.73;
mH = 0.55;
if(d50mm>=0.22)
   fac = (d50mm-0.22)/(0.3-0.22);
   mH = mH+0.45*fac;
   mL = mL+0.27*fac;
   if(d50mm>=0.3)
      mH = 1.;
      mL = 1.;
   end
end
% smooth transition to upper flat bed regime
nH = 1.;
if(Psi>190.)
   nH = 0.5*(1.+cos( pi*(Psi-190.)/(240.-190.) ));
   if(Psi>250.)
      nH = 0.;

   end
end
nL = nH;
Hoa = max(0., mH*nH*(0.275-0.022*Psi.^0.42) );
Loa = max(0., mL*nL*(1.97-0.44*Psi.^0.21) );
return