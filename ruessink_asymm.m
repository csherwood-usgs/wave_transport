function rp = ruessink_asymm( Ur )
% Calculate asymmetry parameters from Ursell number
% rp = ruessink_asymm( Ur )
%
% Input:
%   Ur - Ursell number []
% Returns structure rp:
%   rp.r  - Abreu r (magnitude)
%   rp.phi -Abreau phi (phase) (Eqn. 12)
%   rp.B  - total non-linearity B (Eqn. 7)
%   rp.b  - Malarky & Davies asymmetry (after Eqn 11)
%   rp.Su - Velocity skewness (Eqn. 5)
%   rp.Au - Wave asymmetry (Eqn. 5 with uw replaced with aw)
%
% Ruessink et al., 2012, Coastal Engineering 65:56-63.

% log10 used in Eqn. 9 per G. Ruessink
% Bug in RRvR eqn. 11 solution found by S. Suttles fixed here
% Not sure if next line is right!

dtr = pi/180.;
% Calculate B and phi from RRvR p. 58
p1 = 0.;
p2 = 0.857;
p3 = -0.471;
p4 = 0.297;
B = p1 + (p2 - p1)./(1 + exp( (p3-log10(Ur))./p4 )); % RRvR Eqn. 9
p5 = 0.815;
p6 = 0.672;
psi = dtr*(-90.) + dtr*90. * tanh(p5./(Ur.^p6)); % RRvR Eqn. 10.
% b can be used directly in MD equations
b = sqrt(2.)*B./sqrt(2.*B.^2+9.); % Solution to RRvR Eqn. 11
r = 2.*b./(b.^2+1.);
phi = -psi-pi/2.;                % RRvR Eqn. 12
% dimensionless velocity and acceleration skewness
% RRvR Eqn. 5 and MD Eqn 4a,b
Su = B.*cos(psi); 
Au = B.*sin(psi);

rp.r = r;
rp.B = B;
rp.phi = phi;
rp.b = b;
rp.Su = Su;
rp.Au = Au;
return