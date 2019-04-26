% u3_v_rphi - What is the relationship between asymmetry parameters and
% u^3?

Ur = logspace(-1,1,20);
for i=1:length(Ur)
   rp = ruessink_asymm( Ur(i) )
   
   figure(1); clf
   T = 10
   Uw = 1.
   auw = abreu_uw( rp.r, rp.phi, Uw, 10, 1, 100 );
   u = auw.u;
   u3bar = mean(auw.u.^3);
   urms = rms(u);
   S = u3bar/urms;
   ts = sprintf('Ur=%.4f, r=%.1f, phi=%.1f, ub=%.2f, u^3=%.2f, urms=%.2f, S=%.2f',Ur(i),rp.r,rp.phi,Uw,u3bar,urms,S)
   title(ts)
   pause
end
%%
rlist = [0. : .05 : 1];
nr = length(rlist)
philist =[0.];
np = length(philist);

Uw = 50.
T = 10;

u3 = NaN*ones(np*nr,4);
for i=1:length(rlist)
   for j=1:length(philist)
      r = rlist(i);
      phi = philist(j);
      auw = abreu_uw( r, phi, Uw, T, 0, 100 );
      u3(i,1)=r;
      u3(i,2)=phi;
      u3bar = mean(auw.u.^3)
      u3(i,3)=u3bar;
      u3(i,4)=u3bar/rms(u);
   end
end
%%
figure(1); clf
plot3(u3(:,1),u3(:,2),u3(:,3),'.')
xlabel('r')
ylabel('phi')
zlabel('u^3')