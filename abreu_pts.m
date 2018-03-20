function af = abreu_pts( r, phi, Uw, T )
% Calculate umax, umin, and phases of asymmetrical wave orbital velocity
% af = abreu_pts( r, phi, Uw, T )
% 
% Input:
%   r and phi - asymmetry parameters from Ruessink et al.
%   Uw - amplitude of orbital velocity (m/s)
%   T  - wave period
% 
% Returned:
%    af.T  - wave period (s) 
%    af.DTc  - duration under crest (s)
%    af.DTt  - duration under trough (s)
%       (T = DTc + DTt)
%    af.DTcu - duration of acceleration under crest (s)
%    af.DTtu - duration of acceleration under trough (s)
%    af.Tzd  - time of zero down-crossing (s)
%    af.Tzu  - time of zero up-crossing (s)
%    af.Tc   - time of maximum velocity under crest (s)
%    af.Tt = - time of minimum velocity under trough (s)
%    af.umax - maximum velocity under crest (m/s)
%    af.umin - minimum velocity under trough (m/s)
%    af.R  - Velocity skewness parameter == umax/(umax-umin) ()
%    af.Beta - Acceleration skewness parameter == amax/(amax-amin) ()
%    af.Sk - van Rijn et al. 2011 (Eqn 2) assymetry statistic ()
%    af.As - van Rijn et al. 2011 (Eqn 2) assymetry statistic ()
% 
% Abreu, T., Silva, P.A., Sancho, F., and Temperville (2010) 
%   Analytical approximate wave form for asymmetric waves. 
%   Coastal Engineering, 57(7):656-667.
%   doi: http://dx.doi.org/10.1016/j.coastaleng.2010.02.005.
% 
% Malarkey, J. and A. G. Davies (2012) Free-stream velocity descriptions 
%   under waves with skewness and asymmetry.
%   Coastal Engineering, 68:78-95.
%   doi: http://dx.doi.org/10.1016/j.coastaleng.2012.04.009.
%   
% Ruessink, B. G., G. Ramaekers, and L. C. van Rijn (2012)
%   On the parameterization of the free-stream non-linear wave orbital 
%   motion in nearshore morphodynamic models. Coastal Engineering, 6556-63.
%   doi: http://dx.doi.org/10.1016/j.coastaleng.2012.03.006.
%   
% van Rijn, L. C., P. K. Tonnon, and D. J. R. Walstra (2011)
%   Numerical modelling of erosion and accretion of plane sloping beaches 
%   at different scales. Coastal Engineering, 58: 637-655. 
%   doi: http://dx.doi.org/10.1016/j.coastaleng.2011.01.009.
%
% csherwood@usgs.gov

w = 2*pi/T;

% alternative formulation Eqns 16a,b in Malarkey & Davies
phi = -phi; %
P = sqrt(1.-r*r); % same as f
b = r/(1.+P);

% Appendix E of Malarkey & Davies
% Phase of umax (crest) and umin (trough) (in radians, from 0 to 2*pi)
c = b*sin(phi);
tmc = asin((4.*c*(b*b-c*c)+(1.-b*b)*(1.+b*b-2.*c*c))/((1.+b*b).^2-4.*c*c));
tmt = asin((4.*c*(b*b-c*c)-(1.-b*b)*(1.+b*b-2.*c*c))/((1.+b*b).^2-4.*c*c));
if(tmt<0.)
   tmt = tmt+2.*pi;
end
if(tmc<0.)
   tmc = tmc+2*pi;
end

% umax and umin - non dimensional
umax = 1+c;
umin = umax-2;
% dimensional
umax = Uw*umax;
umin = Uw*umin;

% phase of zero upcrossing and downcrossing (radians)
tzu = asin(b*sin(phi)); % = arcsin(c)
tzd = 2.*acos(c)+tzu;

% Calculate assymetry parameters R and Beta
% R = umax/(umax-umin) % gives same result as:
R = 0.5*(1.+ b*sin(phi) ); % MD Eqn 17
% after MD Eqn. 18
Fo = (r<=.5) * (1.-0.27*(2.*r)^2.1) + ...
     (r> .5) * (0.59 + 0.14*(2.*r)^-6.2);
% MD Eqn. 15b
Br0 =(r< 0.5)* (0.5*(1+r)) + ...
     (r>=0.5)* ( 4.*r*(1.+r)/(4.*r*(1.+r)+1.) );
 
% next eqn. corrected to match Malarkey/Davies(2012) Eqn 18 CRS 3/19/2018 
Beta = 0.5+(Br0-0.5)*sin((0.5*pi-abs(phi))*Fo)/sin(0.5*pi*Fo);
  
% Calculate assymetry parameters Sk and As (same as van Rijn et al. 2011 Eqn. 2)
% next two eqns. corrected to match Malarkey/Davis(2012) Table 1. CRS
% 3/19/2018
Sk =  3.*b*sin(phi)/sqrt(2.*(1.+b^2)^3);
As = -3.*b*cos(phi)/sqrt(2.*(1.+b^2)^3);

% Could also use MD Appendix C to calculate uspike, aspike, and other
% measures of assymetry

% These are the dimensional fractions of wave periods needed by Van der A eqn.
af.T = T;
af.DTc  = (tzd-tzu)/w;
af.DTt  = T - af.DTc;
af.DTcu = (tmc-tzu)/w;
af.DTtu = (tmt-tzd)/w;

af.Tzd = (tzd)/w;
af.Tzu = (tzu)/w;
af.Tc =  (tmc)/w;
af.Tt  = (tmt)/w;

af.umax = umax;
af.umin = umin;
af.R = R;
af.Beta = Beta;
af.Sk = Sk;
af.As = As;

return