function wtout = wave_transport_func( wtin )
% Compute wave transport under asymmetrical waves
%
% References:
% RRvR = Ruessink, Ramaekers, and van Rijn (2012), Coastal Eng. 65:56-63
% A13  = van der A et al., (2013), Coastal Eng. 76:26-42.
%
% Calls:
%   fw_func
%   fd_func
%   dsf_func
%   ksw_func
%   ksd_func
%   od_ripple
%
% csherwood@usgs.gov
% 13 July 2015
% Input:
%   Hs - sig. wave height (m)
%   Td - dom. wave period (s)
%    h - water depth (m)
Hs = wtin.Hs;
Td = wtin.Td;
h =  wtin.h;
if(length(wtin.d)==1)
d50 = wtin.d;
d90 = 1.5*d50; % check this...van Rijn has an eqn.
else
   % TODO - treat multiple size classes
   error('single grain size only')
end
% constants
dtr = pi/180.; % degrees to radians
g = 9.81;
vk = 0.41;
rhow = 1027.;
rhos = 2650.;
nu = 1.36e-6;

% Sed properties - need tau_crit and ws
% s = (rhos-rhow)/rhow after A13, Eqn. 1, but should be s = rhos/rhow
% (Soulsby, 1997)
s = rhos/rhow;
p = soulsby_particle(d50,rhos,rhow,nu)
theta_crit = p.Shields_crit;
tau_crit = theta_crit * (g*(s-1.)*d50) % Soulsby Eqn 74 (inverse)
ws = p.ws;

% TODO - Should this be calculated from Hs?
[uhat,Tbav]=ubspecfun(Hs,Td,h);
%uhat = 0.5;    % "significant" orbital velocity (ubr = ubot, or sqrt(2) times that?
ahat = uhat*Td/(2.*pi);

% Current speed and direction
mag_u_d = wtin.mag_u_d;
dir_u_d = wtin.dir_u_d;  % direction of current, meas. CCW from wave direction (A13 Fig. 2)
costhet = cos(pi*dir_u_d/180.);
sinthet = sin(pi*dir_u_d/180.);

% Max. mobility number (Appdx. B) for irregular waves...
% but want to use max(crest, trough) later
Psi = (1.27*uhat)^2. / (g*(s-1.)*d50);
[Hoa, Loa] = od_ripple( d50, Psi );
rh = Hoa*ahat
rl = Loa*ahat

% Calculate wave-averaged stress and roughness.
% Iteration is required because roughness varies with stress
dsf = d50; % first guess: no sheet flow layer
ksw = ksw_func( d50, rh, rl, 0. );
ksd = ksd_func( d50, d90, rh, rl, 0. );

tlast = 999.;
acc = 99.;
tol = 1e-3;
nit = 20;
i = 0;
sf = 1./((s-1.)*g*d50);
while (i<nit && (abs(acc)>tol) )
   fw = fw_func( ahat, ksw );
   fd = fd_func(  dsf, ksd );
   mean_mag_theta = sf*(0.5*fd*mag_u_d^2. + 0.25*fw*uhat^2. );
   acc = 2.*(tlast-mean_mag_theta)/(eps+tlast+mean_mag_theta);
   tlast = mean_mag_theta;
   fprintf(1,'%02d fd=%6.4f fw=%6.4f theta=%6.3f acc =%6.4f dsf=%6.4f ksw=%6.4f ksd=%6.4f\n',...
      i, fd, fw, mean_mag_theta, acc, dsf, ksw, ksd);
   dsf = dsf_func( d50, mean_mag_theta );
   ksw = ksw_func( d50, rh, rl, mean_mag_theta );
   ksd = ksd_func( d50, d90, rh, rl, mean_mag_theta );
   i=i+1;
end

% wavenumber k
w = 2.*pi/Td;    % angular wave frequency
kh = qkhfs(w,h); %
k = kh/h;        % wave number
cw = 2*pi/(k*Td); % wave celerity
fprintf(1,'Wave properties: \n Td: %6.2f s\n w : %6.4f\n k : %6.4f\n cw: %6.2f m/s\n',...
   Td,w,k,cw)

% Ursell number Ur
Ur = 0.75*0.5*Hs*k./(kh.^3) % RRvR Eqn. 6.
fprintf(1,'Ruessink asymmetery parameters: \n')
rp = ruessink_asymm( Ur )
r = rp.r;
phi = rp.phi;

% Or, specify to compare w/ Abreu figures
% r = .8; phi=0.

% Abreu orbital velocity time series can be computed from r and phi
% according to RRvR Eqn. 4, which is same as Abreu Eqn. 7.
% This routine provides the values you need for the transport calcs:
fprintf(1,'Points in velocity time series from Abreu: \n')
af = abreu_pts(r, phi, uhat, Td );
% key intervals during wave period...same notation at A13, Fig 1.
T = af.T
Tc = af.DTc
Tcu = af.DTcu
Tt = af.DTt
Ttu = af.DTtu
R = af.R

% This is a time series...not necessary, but use to check:
if(0)
   figure(1); clf; hold on
   auw = abreu_uw( r, phi, Uw, Td, 1, 101 )
end

% Crest and trough "representative" velocities:
% (coordinate system is aligned with wave direction; A13, Fig. 2)
% TODO - check to make sure these orbital velocities are correct and don't
% have to be made into "representative" by multiplying by, say, sqrt(2)
uhatc =  af.umax;  % see text after A13 Eqn. 7
uhatt = -af.umin;  % I think minus sign is needed to offset minus sign in Eqn. 13 below
utildecr = 0.5*sqrt(2)*uhatc; % A13 Eqn. 10
utildetr = 0.5*sqrt(2)*uhatt; % A13 Eqn. 11

% Crest and trough velocity components and magnitude^2:
ucrx = utildecr+mag_u_d*costhet;    % A13 Eqn. 12
ucry = mag_u_d*sinthet;
mag_ucr = sqrt( ucrx^2 + ucry^2 );
utrx = -utildetr + mag_u_d*costhet; % A13 Eqn. 13
utry = mag_u_d*sinthet;
mag_utr = sqrt( utrx^2+utry^2 );
fprintf(1,'Crest and trough velocity components: \n')
fprintf(1,'  ucrx, ucry: %f %f\n', ucrx, ucry)
fprintf(1,'  utrx, utry: %f %f\n', utrx, utry)

% Check to see if R and Beta calculated from these matches
% A13, text after Eqn. 13)
Rcheck = uhatc/(uhatc+uhatt)

%% Calculate Shields parameters
% This requires calculating magnitude for crest and trough according
% to A13 Eqns. 13 - 21.
alpha = mag_u_d/(mag_u_d + uhat); %A13 Eqn. 19
fprintf(1,'Current fraction: %f\n',alpha);
sf = 1./((s-1)*g*d50);

% Wave-averaged Reynolds stress from streaming (A13 Eqn. 22)
% (does not occur in u-tube experiments)
alphaw = 4. /(3.*pi);
fwd = alpha*fd + (1.-alpha)*fw;
tauwre = rhow*fwd/(2.*cw)*alphaw*uhat^3

c1 = 2.6; % coefficient in friction factor; A13, text after Eqn 21

% Crest
fwc = 0.3;
if( ahat / ksw ) > 1.587
   fwc = 0.00251*exp( 5.21* ((ahat*(2*Tcu/Tc)^c1)/ksw)^(-0.19) );
end
fwdc = alpha*fd + (1-alpha)*fwc; %A13 Eqn. 18
Sc = sf*(0.5*fd*mag_u_d^2. + 0.25*fwdc*mag_ucr^2 );
Scx = Sc*ucrx/mag_ucr + sf*tauwre/rhow; %A13 Eqn. 15 (see note below)
Scy = Sc*ucry/mag_ucr;             %A13 Eqn. 16
% Note on the rhow...this seems to be required to get dimensions correct; s
% has no units, but Tau is in Pascals. Related to (I think) mistaken
% definition of s after A13, Eqn 1.

% Trough
fwt = 0.3;
if( ahat / ksw ) > 1.587
   fwt = 0.00251*exp( 5.21* ((ahat*(2*Ttu/Tt)^c1)/ksw)^(-0.19) );
end
fwdt = alpha*fd + (1-alpha)*fwt; %A13 Eqn. 18
St = sf*(0.5*fd*mag_u_d^2. + 0.25*fwdt*mag_utr^2 );
Stx = St*utrx/mag_utr + sf*tauwre/rhow; %A13 Eqn. 15 (see not above)
Sty = St*utry/mag_utr;             %A13 Eqn. 16

fprintf(1,'For comparision wiht A13 Fig. 3: \n')
fprintf(1,'  ahat/ksw : %6.0f \n   Beta   : %4.1f\n',ahat/ksw,af.Beta)
fprintf(1,'  Crest:  fwdc and Shields param: %f, %f\n',fwdc,Sc)
fprintf(1,'  Trough: fwdt and Shields param: %f, %f\n',fwdt,St)
fprintf(1,'thetac / thetat : %4.2f\n',Sc/St)


%% adjusted settling velocities
% assume suspended sediment is smaller than d50
ds = 0.8*d50;
padj = soulsby_particle(ds,rhos,rhow,nu);
wsa = padj.ws;
% correct for vertical flow using second-order Stokes waves at elevation of
% ripple crests or sheet-flow layer thickness
zws = max(rh,dsf);
% Following calcs of vertical velocity under crest and trough
% are poached from santoss_core.m but seems to underestimate for really
% asymmetrical waves, compared with mins and maxes from time series\
% I calculate in stokes_second_order.m. Same for celerity: it matches my
% calcs until starting to push H/h
% I am not sure why both crest and trough are calculated...they are
% exactly the same. Not sure they should be, but that is also what I get in
% stokes_second_order.m
worbc1= pi*Hs*zws/Td/h;
worbt1= pi*Hs*zws/Td/h;
worbc2= worbc1*2*(R+R-1);
worbt2= worbt1*2*(R+R-1);
worbc = (1/8)*worbc1*sqrt(64-(-worbc1+sqrt(worbc1^2+32*worbc2^2))^2/(worbc2^2))+...
   worbc2*sin(2*acos((1/8)* (-worbc1+sqrt(worbc1^2+32*worbc2^2))/worbc2));
worbt = (1/8)*worbt1*sqrt(64-(-worbt1+sqrt(worbt1^2+32*worbt2^2))^2/(worbt2^2))+...
   worbt2*sin(2*acos((1/8)* (-worbt1+sqrt(worbt1^2+32*worbt2^2))/worbt2));
wsc = max(0.,wsa+worbc) % enhanced settling
wst = max(0.,wsa-worbt) % inhibited settling
%%
% wave celerity for second-order Stokes wave, Table 5.3 in Komar, 1998.
L = 2*pi/k;
term2 = (pi*Hs/L)^2 *(5. + 2.*cosh(4*pi*h/L) + 2.*(cosh(4.*pi*h/L))^2) /(8.*(sinh(2*pi*h/L))^4);
C = (g/w)*tanh(kh)*(1.+ term2); % celerity
xi = 1.7;     % I think...from santoss_core.m
alphar = 8.2; % calibration term alpha - A13, Section 2.5
Pc = alphar*(1-xi*uhatc/C)*zws/(wsc*2*Tcu);
Pt = alphar*(1+xi*uhatt/C)*zws/(wst*2*Ttu);
%% Compute flux magnitudes
% m and n are key calibration parameters...A13, bottom of section 2.5
m = 11.0;
n = 1.2;
Oc = m*max(Sc-theta_crit,0.)^n; % A13, Eqn 2.
if(Pc<=1.) % A13 Eqn 23, 24
   Occ = Oc;
   Oct = 0.;
else
   Occ = Oc/Pc;
   Oct = (1.-1./Pc)*Oc;
end
Ot = m*max(St-theta_crit,0.)^n; % A13, Eqn 2.
if(Pt<=1.)  % A13 Eqn 23, 26
   Ott = Ot;
   Otc = 0.;
else
   Ott = Ot/Pt;
   Otc = (1.-1./Pt)*Ot;
end
% Non-dimensional fluxes from A13, Eqn. 1
Phicx = sqrt(Sc)*Tc/T*(Occ+Tc/(2*Tcu)*Otc)*Scx/Sc
Phicy = sqrt(Sc)*Tc/T*(Occ+Tc/(2*Tcu)*Otc)*Scy/Sc
Phitx = sqrt(St)*Tt/T*(Ott+Tt/(2*Ttu)*Oct)*Stx/St
Phity = sqrt(St)*Tt/T*(Ott+Tt/(2*Ttu)*Oct)*Sty/St
% Multiply by sfd3 to get dimensional fluxes
sfd3 = sqrt((s-1)*g*d50^3);
wtout.Ur = Ur;
wtout.R = af.R;
wtout.Beta = af.Beta;
wtout.fx = sfd3*(Phicx+Phitx);
wtout.fy = sfd3*(Phicy+Phity);
