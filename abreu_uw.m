function auw = abreu_uw( r, phi, Uw, T, iplot, n )
% Time series of assumetric wave-orbital velocties and accelerations
% auw = abreu_uw( r, phi, Uw, T, iplot, n )
%
% Input:
%   r and phi - assymetry parameters from Ruessink et al. 2012
%   Uw - amplitude of orbital velocity (m/s)
%   T  - wave period (s)
%   iplot - 0 or 1: option to make plot (default = 0 )
%   n - number of points in time series (default = 50)
iplot = 1;
if(exist('iplot','var')~=1),iplot=0;,end
if(exist('n','var')~=1),n=50;,end;

w = 2.*pi/T;
wt = linspace(0,2.*pi,n)';
% Abreu eqn., also MD12 eqns. 13 a,b
f = sqrt( 1. - r^2 )
numer = sin(wt) + ( r*sin(phi) ./ (1.+sqrt(1.-r.^2)) );
denom = (1.-r*cos(wt+phi));
ut = Uw*f*numer ./ denom;
numer2 = cos(wt)-r*cos(phi)-r.^2 ./ (1.+sqrt(1.-r.^2))*sin(phi)*sin(wt+phi);
at = Uw*w*f*numer2 ./ denom.^2;
t = wt/w;
u = ut/Uw;
a = at/(Uw*w;
if(iplot)
   % plot time series of velocity and acceleration
   h2=plot(t,a,'-r','linewidth',2);
   hold on
   h1=plot(t,u,'-b','linewidth',2);
   ylabel('u (m/s); du/dt (m/s^2)')
   xlabel('Time (s)')
   % find critical points and plot
   af = abreu_pts( r, phi, Uw, T );
   plot(af.Tc,af.umax,'ob')
   plot(af.Tt,af.umin,'ob')
   plot(af.Tt,0,'+k')
   plot(af.Tc,0,'+k')
   plot(af.Tzd,0,'*k')
   plot(af.Tzu,0,'*k')
   h3=plot([af.Tzu, af.Tzu+af.DTc], [-1.2 -1.2],'-c','linewidth',5);
   h4=plot([af.Tzd, af.Tzd+af.DTt], [-1.2 -1.2],'-m','linewidth',5);
   h5=plot([af.Tzu, af.Tzu+af.DTcu], [-1.1 -1.1],'--c','linewidth',5);
   h6=plot([af.Tzd, af.Tzd+af.DTtu], [-1.1 -1.1],'--m','linewidth',5);
   legend([h1;h2;h3;h4;h5;h6],'Velocity','Acceleration','T_c','T_t','T_{cu}','T_{tu}')
   grid on
end
auw.t = t;
auw.u = u;
auw.a = a;
return



