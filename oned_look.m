
a=load('./swan_shoal/CRS1d.dat');
x = a(:,1);
h = a(:,2);
h(h<=-99.00)=NaN;
Hs = a(:,3);
Hs(Hs<=-9.00)=NaN;
setup = a(:,10);
setup(setup<=-9.00)=NaN;
Tm01 = a(:,4);
Tm01(Tm01<=-9.00)=NaN;
Ubot = a(:,11);
Ubot(Ubot<=-9.00)=NaN;
Tmbot = a(:,12);

slope = abs(diff(h(1:2))/diff(x(1:2)))
g = 9.81;
Lo = Tm01(1)^2 *g/(2*pi)
Ib = slope/sqrt(Hs(1)/Lo)

k = NaN*ones(size(Tm01));
Ur = NaN*ones(size(Tm01));
r = NaN*ones(size(Tm01));
phi = NaN*ones(size(Tm01));

for i=1:length(k)
    k(i)=qkhfs( 2*pi/Tm01(i), h(i))/h(i);
    Ur(i) = 0.75*0.5*Hs(i)*k(i)./(k(i)*h(i).^3); % RRvR Eqn. 6.
    rp = ruessink_asymm( Ur(i) );
    r(i) = rp.r;
    phi(i) = rp.phi;
end
%%
figure(1)
clf
subplot(411)
plot(x,zeros(size(x)),'--k')
hold on
plot(x,-h,'-k','linewidth',3)
hold on
plot(x,Hs,'-b','linewidth',2)
plot(x,setup,'-m','linewidth',2)

subplot(412)
plot(x,r)
hold on
plot(x,phi)

subplot(413)
plot(x,Ubot)
