% run_wave_transport - Script to run wave_transport_func
% 
wtin.Hs = 2;
wtin.Td = 10.;
wtin.d = 0.15e-3;
wtin.mag_u_d = 0.5
wtin.dir_u_d = 90.  % direction of current, meas. CCW from wave direction (A13 Fig. 2)

slope = 1/20
ho = -5*wtin.Hs;
xo = ho/slope;
x = linspace(xo,Hs/0.8,100);
h = x*slope;
%%
for i=1:length(h)
   wtin.h = -h(i);
   wtout(i)=wave_transport_func(wtin)
end

