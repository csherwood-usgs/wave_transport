
% P. D. Komar, "Beach Processes and Sedimentation", Second Edition.
% Table 5.3 (with corrections)
H = 1;
T = 12;
h = 5;
g = 9.81;
zup = 0.1;    % elevation of u and w above seafloor
z = -(h-zup);
w = 2*pi/T;
kh = qkhfs(w,h);
k = kh/h;
L = 2*pi/k;
nx = 9;
nt = 49;
x = L*(0:nx)'/nx;
t = T*(0:nt)/nt;
x = repmat(x,1,nt+1);
t = repmat(t,nx+1,1);
% check validity per Fig. 5-29 in Komar
if(H/h>0.78),warning('H/h>0.78'),end
if(H/L>0.142*tanh(kh)),warning('H/L>0.142*tanh(kh)'),end
%% From Svendsen
kxwt = k*x+w*t;
eta = 0.5*H*cos(kxwt) + (1./16.)*k*H.^2*(3*coth(kh).^3-coth(kh))*cos(2*kxwt); % (8.2.76)
%% From Komar
eta_komar = 0.5*H*cos(k*x+w*t) + 0.5*pi*H*H/L*cosh(kh)*(2+cosh(2*kh))/(4*sinh(kh)^3)*(cos(2*(k*x-w*t))); %(Table 5.3, corrected)
%%
figure(1); clf
subplot(311)
plot(t(1,:),eta(1,:))
hold on
plot(t(1,:),eta_komar(1,:),'.k')
ylabel('Height (m)')
shg
%%
term2 = (pi*H/L)^2 *(5. + 2.*cosh(4*pi*h/L) + 2.*(cosh(4.*pi*h/L))^2) /(8.*(sinh(2*pi*h/L))^4);
C = (g/w)*tanh(kh)*(1.+ term2); % celerity
au = ((pi*H/T)*cosh(k*(z+h))/sinh(kh))
bu = ((3./4.)*(pi*H/L)^2*C*cosh(2*k*(z+h))/sinh(kh)^4)
u = au*cos(k*x+w*t) + bu*cos(2*(k*x-w*t));
umax = max(u(1,:))
umin = min(u(1,:))
R = umax/(umax-umin)

aw = ((pi*H/T)*sinh(k*(z+h))/sinh(kh))
bw = ((3./4.)*(pi*H/L)^2*C*sinh(2*k*(z+h))/sinh(kh)^4)
wwa = aw*sin(k*x-w*t);
wwb = bw*sin(2*(k*x-w*t));
ww = aw*sin(k*x-w*t)+ bw*sin(2*(k*x-w*t));

fprintf('Min and max w: %6.4f %6.4f\n',min(ww(1,:)),max(ww(1,:)))

figure(1);
subplot(312)
h2=plot(t(1,:),u(1,:));
ylabel('u (m/s)')
subplot(313)
h3=plot(t(1,:),ww(1,:));
ylabel('v (m/s)')
hold on
h4 = plot(t(1,:),wwa(1,:));
h5 = plot(t(1,:),wwb(1,:));
xlabel('t (s)')
%legend([h1;h2;h3],'eta','u','w')
legend([h3;h4;h5],'w','aw','bw')
ts = sprintf('H = %4.1f, T = %4.1f, h = %5.1f, zup = %4.2f',H,T,h,zup);
title(ts)

sq = 4*bu^2-au^2
if(sq)<0,sq=0;end
tu = atan2(au/(2*bu),-sqrt(sq)/(2*bu))

% Following calcs of vertical velocity under crest and trough
% are poached from santoss_core.m but seems to underestimate for really
% asymmetrical waves, compared with mins and maxes from time series\
% I calculate in stokes_second_order.m. Same for celerity: it matches my
% calcs until starting to push H/h 

worbc1= pi*H*zup/T/h;
worbt1= pi*H*zup/T/h;
worbc2= worbc1*2*(R+R-1);
worbt2= worbt1*2*(R+R-1);

worbc = (1/8)*worbc1*sqrt(64-(-worbc1+sqrt(worbc1^2+32*worbc2^2))^2/(worbc2^2))+...
   worbc2*sin(2*acos((1/8)* (-worbc1+sqrt(worbc1^2+32*worbc2^2))/worbc2))

worbt = (1/8)*worbt1*sqrt(64-(-worbt1+sqrt(worbt1^2+32*worbt2^2))^2/(worbt2^2))+...
   worbt2*sin(2*acos((1/8)* (-worbt1+sqrt(worbt1^2+32*worbt2^2))/worbt2))

% Calculation wave propagation speed c
KSI = 4*pi^2*h/(g*T^2)  % (h=depth, T=wave period)
   e = sqrt(KSI) * (1. + 0.2 * KSI);
if KSI > 1.
   e = KSI * (1.+0.2*exp(2. - 2*KSI));
end
L2= 2* pi* h / e;
c = L2 / T;
fprintf('CRS celerity: %f6.3; SANTOSS celerity: %6.3f\n',C,c);
