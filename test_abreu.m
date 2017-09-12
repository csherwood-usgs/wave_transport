% test_abreu - Check to see the Abreu time series works
% select r and phi to match curves in Figs 1 and 2
philist = [0. -pi/4 -pi/2]
rlist = (0.:.25:.75)
for j=1:length(philist)
   for i=1:length(rlist)
      
      phi = philist(j)
      r = rlist(i)
      
      Uw = 50.
      T = 14.
      
      af = abreu_pts(r, phi, Uw, T)
      
      
      n = 50
      w = 2.*pi/T;
      wt = linspace( 0., 2.*pi, n); % phase
      % Abreu eqn., also MD12 eqns. 13 a,b
      f = sqrt( 1. - r.^2 );
      numer = sin(wt) + ( r*sin(phi)/(1.+sqrt(1.-r.^2)) );
      denom = (1.-r*cos(wt+phi));
      ut = Uw*f*numer./denom;
      numer2 = cos(wt)-r*cos(phi)-r.^2/(1.+sqrt(1.-r.^2))*sin(phi)*sin(wt+phi);
      at = Uw*w*f*numer2./denom.^2;
      
      % alternative formulation Eqns 16a,b in Malarkey & Davies
      phi = -phi; %
      P = sqrt(1.-r*r); % same as f
      b = r/(1.+P)
      
      fbt = 1.-b*b;
      numer = sin(wt)-b*sin(phi);
      denom = (1.+b*b-2.*b*cos(wt-phi));
      utm = Uw*fbt*numer./denom;
      numer2 = (1.+b*b)*cos(wt)-2.*b*cos(phi)+2.*b*b*sin(phi)*sin(wt-phi);
      atm = Uw*w*fbt*numer2./denom.^2;
      
      % Appendix E of Malarkey * Davies
      % Phase of umax, and umin
      c = b*sin(phi)
      tmm = asin((4.*c*(b*b-c*c)-(1.-b*b)*(1.+b*b-2.*c*c))/((1.+b*b).^2-4.*c*c));
      tmp = asin((4.*c*(b*b-c*c)+(1.-b*b)*(1.+b*b-2.*c*c))/((1.+b*b).^2-4.*c*c));
      if(tmm<0.)
         tmm = tmm+2.*pi;
      end
      if(tmp<0.)
         tmp = tmp+2*pi;
      end
      % non dimensional
      umax = 1+c;
      umin = umax-2;
      % dimensional
      umax = Uw*umax
      umin = Uw*umin
      
      % zero upcrossing
      tz = asin(b*sin(phi)); % = arcsin(c)
      tzd = 2.*acos(c)+tz;
      % sigma Eqn 19
      sig1 = asin( (4.*c*(b*b-c*c)+(1.-b*b)*(1.+b*b-2.*c*c))/((1.+b*b).^2-4.*c*c) );
      if( phi <= 0.5*pi )
         sig = (1./pi)*(sig1-tz);
      else
         sig = (1./pi)*(pi-sig1-tz);
      end
      (tmp-tz)/pi; % sigma from Eqn 5
      
      %These are the dimensional fractions of wave periods needed by Van der A eqn.
      Tcu = (tmp-tz)/w
      Tc  = (tzd-tz)/w
      Ttu = (tmm-tz)/w
      Tzd = (tzd-tz)/w
      Tt = T-Tc
      
      figure(j); if(i==1),clf,end;
      plot(wt/(2.*pi),ut/Uw)
      hold on
      plot(wt/(2.*pi),at/(Uw*w))
      plot(wt/(2.*pi),utm/Uw,'--')
      plot(wt/(2.*pi),atm/(Uw*w),'--')
      plot(tmm/(2.*pi),umin/Uw,'or')
      plot(tmp/(2.*pi),umax/Uw,'ob')
      plot(tz/(2.*pi),0,'ok')
      plot(tzd/(2.*pi),0,'ok')
      xlabel('{\itt/T}')
      ylabel('{\itu/U_w,  a/(U_w \omega)}')
   end
end