function [R,Beta] = R_Beta_from_r_phi(r,phi)
% [R,Beta] = R_Beta_from_r_phi(r,phi)
%  Recover velocity asymmetry and accel. assymetry from Reussink params
%

P = sqrt(1-r^2);
b = r/(1.+P);
% Malarky & Davies Eqn 17
R = 0.5*(1+b*sin(phi));

if (r >= 0) && (r <= 0.5)
    Betar0 = 0.5*(1.+r);
elseif ((r >= 0.5) && (r < 1))
    Betar0 = 4.*r*(1+r)/(4.*r*(1+r)+1.);
else
    error('r out of range')
    Beta = NaN;
    return
end
if(r<=0.5)
    F0 = 1.-0.27*(2.*r)^2.1;
elseif r > 0.5
    F0 = 0.58 + 0.14*(2.*r)^(-6.2);
end
Beta = 0.5 + (Betar0 - 0.5)*sin(0.5*pi-abs(phi)*F0)/sin(0.5*pi*F0);
return
end

