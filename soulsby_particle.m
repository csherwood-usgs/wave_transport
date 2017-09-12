function p = soulsby_particle( d, rhos, rhow, nu, iverbose )
% SOULSBY_PARTICLE - Critical shear stress and settling velocity (MKS units)
% p = soulsby_particle( d, rhos, rhow, nu, iverbose )
%
% Input (MKS units):
%   d - diameter (m)
%   rhos - sediment density (optional, default = 2650 kg/m3)
%   rhow - water density (optional, default = 1027 kg/m3)
%   nu   - kinematic viscosity (optional, default = 1.36e-6 m2 s-1)
%
% From Soulsy (1997) 'Dynamics of Marine Sands'

% csherwood@usgs.gov
% last revised February 2010
if(exist('rhos')~=1),
    rhos = 2650; % particle density (quartz) km m-3
end
if(exist('rhow')~=1),
    rhow = 1027; % water density kg m-3
end
if(exist('nu')~=1),
    nu = 1.36e-6; % kinematic viscosity m2 s-1
end
if(exist('iverbose')~=1),
   iverbose = 0;
end

g = 9.81; % m s-2
s = rhos/rhow; % p. 104
Dstar = d.*(g*(s-1)./(nu.^2)).^(1/3); % Eqn 75
Shields_crit = 0.3./(1+1.2*Dstar) + 0.055*(1-exp(-0.020*Dstar)); % Eqn 77
tau_crit = Shields_crit .* (g*(rhos-rhow).*d ); % Eqn 74 (inverse)
D3 = Dstar.^3;
ws_h = nu*D3                ./(18*d) .* (D3 <= 39) + ...
        nu*Dstar.^(2.1)      ./(6*d)  .* ((D3 > 39) & (D3 < 1e4 )) + ...
        nu*1.05*Dstar.^(1.5) ./d      .* ((D3 >= 1e4) & (D3 < 3e6 )); % Eqn 100
ws_vr = nu*D3                ./(18*d) .* (D3 <= 16.187) + ...
        (nu*10 ./d) .*(sqrt(1+0.01*D3)-1).* ((D3 > 16.187) & (D3 <= 16187 )) + ...
        nu*1.1*Dstar.^(1.5) ./d      .* ((D3 > 16187) & (D3 < 3e6 )); % Eqn 101
ws = (nu ./ d ) .* ( sqrt( 10.36^2 + 1.049 * D3) - 10.36 ); % Eqn 102
if(iverbose)
   fprintf(1,'Grain size (phi)                   d = %g [phi]\n',-1.44296504*log(d*1e3));
   fprintf(1,'Grain size (microns)               d = %g [um]\n',d*1e6)
   fprintf(1,'Particle density                rhos = %g [kg m-1]\n',rhos)
   fprintf(1,'Water density                   rhow = %g [kg m-1]\n',rhow)
   fprintf(1,'Kinematic viscosity               nu = %g [m2 s-1]\n',nu)
   fprintf(1,'Dimensionless grain size          D* = %g [ ]\n',Dstar)
   fprintf(1,'Threshold Shields param Shields_crit = %g [ ]\n',Shields_crit)
   fprintf(1,'Critical shear stress       tau_crit = %g [N m-2]\n',tau_crit)
   fprintf(1,'Hallermeier settling velocity     ws_vr = %g [m s-1]\n',ws_h);
   fprintf(1,'van Rijn settling velocity        ws_vr = %g [m s-1]\n',ws_vr);
   fprintf(1,'Soulsby settling velocity            ws = %g [m s-1]\n',ws);
end
p.Dstar = Dstar;
p.Shields_crit = Shields_crit;
p.tau_crit = tau_crit;
p.ws_vr = ws_vr;
p.ws = ws;