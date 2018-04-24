% vandera.m - Script to check asymmetric wave calcs
% This is based on vandera.F, the prototype code intended for ROMS, as of
% March 21, 2018

% MODULE mod_scalars
% note that these are not common to functions, so make sure same values are
% used there
rho0 = 1025.0;
g = 9.81 ;            % m/s2
vonKar = 0.41 ;       % non-dimensional
nu = 1.36E-6  ;       % kinematic viscosity m2 s-1
% pi = 3.14159265358979323846_r8 % use Matlab constant
deg2rad = pi / 180.0;
rad2deg = 180.0 / pi;
% END MODULE mod_scalars
% This coefficient modulates the progressive wave stress tauwRe
% Santoss includes this, and discusses values of 0.55 or zero.
% Using 1 matches equations in VA2013.
eps_eff = 1.
%
%      Input wave parameters
%      Hs=significant wave height
%      Td=dominant wave period
%      depth=water depth(m)
%
Hs = 2.0
Td = 10.0
depth =  10.0
d50 = 0.15e-3
d90 = 1.5*d50
rhos= 2650.0
%
% umag_curr is the current velocity magnitude (MAKE Sure its the magnitude)
% direction of currents is counter-clockwise from wave direction
%
umag_curr=abs(0.2)
%phi_curwave=45.0*deg2rad
phi_curwave=0.0*deg2rad


%
% uhat and ahat for the entire wave cycle  uhat is the wave orbital velocity defiend for the entire wave cycle
% TSK --. HARDWIRED for now but later needs to be based on a function for JONSWAP spectrum
uhat=0.5472
%
ahat=uhat*Td/(2.0*pi)
k=kh_calc(Td,depth)/depth     % Wave number
% compare w/ CRS function
kh_CRS = qkhfs( 2*pi/Td, depth);
k_CRS = kh_CRS/depth;
%
c_w=2*pi/(k*Td) ;         % Wave speed

% VA-2013 equation 1 is solved in 3 sub-steps
%
%-----------------------------------------------------------------------
% Ruessink et al. provides equations for calculating skewness parameters
% Uses Malarkey and Davies equations to get "bb" and "r"
% COMMON TO both crest and trough
%-----------------------------------------------------------------------
%
[r, phi ] = skewness_params(Hs, Td, depth )
% Ursell number Ur
Ur = 0.75*0.5*Hs*k_CRS./(kh_CRS.^3) % RRvR Eqn. 6.
% compare with CRS function
rp_CRS = ruessink_asymm( Ur ); % This agrees
%
%-----------------------------------------------------------------------
% Abreu et al. use skewness params to get representative critical orbital
% velocity for crest and trough cycles , get r and phi from above
%-----------------------------------------------------------------------
%
[T_c, T_t, T_cu, T_tu, umax, umin, RR] = ...
    abreu_points(r, phi, uhat, Td)
af_CRS = abreu_pts(r, phi, uhat, Td) % This agrees
%
%-----------------------------------------------------------------------
%           Crest half cycle
%-----------------------------------------------------------------------
% Get the "representative crest half cycle water particle velocity
%    as well as full cycle orbital velocity and excursion
%-----------------------------------------------------------------------
%
% from Abreu points
uhat_c=umax;
% CRS - we flipped this sign
% uhat_t=-umin;
% I think that is wrong...we want uhat_t to be negative
uhat_t=umin;

%
%-----------------------------------------------------------------------
% VA2013 Equation 10, 11
%-----------------------------------------------------------------------
%
uc_r=0.5*sqrt(2.0)*uhat_c;
ut_r=0.5*sqrt(2.0)*uhat_t;
% VA2013 Eqn 12
uc_rx=uc_r+umag_curr*cos(phi_curwave);
uc_ry=umag_curr*sin(phi_curwave);
mag_uc=sqrt(uc_rx*uc_rx+uc_ry*uc_ry);
fprintf(1,'crest vel. %f, %f\n',uc_rx, uc_ry)
%            alpha=umag_wave/(umag_wave+uhat)
%
%-----------------------------------------------------------------------
%           Trough half cycle
%-----------------------------------------------------------------------
% 1. Get the representative trough half cycle water particle velocity
%    as well as full cycle orbital velocity and excursion
%-----------------------------------------------------------------------
%
% VA2013 Eqn 13
% TODO CRS - this had a minus sign to make ut_ry negative
% ut_rx=-ut_r+umag_curr*cos(phi_curwave)
ut_rx=ut_r+umag_curr*cos(phi_curwave)
ut_ry=umag_curr*sin(phi_curwave)
mag_ut=sqrt(ut_rx*ut_rx+ut_ry*ut_ry)
fprintf(1,'trough vel. %f, %f\n', ut_rx, ut_ry)
%
% This should agree with af_CRS.R 
Rcheck = uhat_c/(uhat_c+uhat_t)
% af_CRS.R % This agrees

%
%-----------------------------------------------------------------------
% 2. Bed shear stress (Shields parameter) for Crest half cycle
%    alpha VA2013 Eqn. 19
%-----------------------------------------------------------------------
% TSK    smgd=(Srho(ised,ng)/rho0-1.0)*g*sd50(ised,ng)
%
smgd=(rhos/rho0-1.0)*g*d50;
osmgd=1.0/smgd;
%
% TSK- Think about are multiple sediments
%
[eta_c, dsf_c, mag_theta_c, tau_wRe_c] = ...
    stress_progressive_surface_waves(d50, d90, osmgd,...
    Td, depth,...
    umag_curr, uhat, uhat_c, ahat, ...
    T_cu, T_c, mag_uc);
%
%-----------------------------------------------------------------------
% Shields parameter in crest cycle VA2013 Eqn 15
% rho0 is required for non-dimensionalizing
%-----------------------------------------------------------------------
%
theta_cx=abs(mag_theta_c)*uc_rx/(mag_uc)+eps_eff*tau_wRe_c*osmgd/rho0;
theta_cy=abs(mag_theta_c)*uc_ry/(mag_uc);
mag_theta_c=sqrt(theta_cx*theta_cx+theta_cy*theta_cy);
fprintf(1,'|theta_c|,theta_cx, theta_cy, theta_cprog: %f %f %f %f\n',...
    mag_theta_c, theta_cx, theta_cy, eps_eff*tau_wRe_c*osmgd/rho0 )

%
%-----------------------------------------------------------------------
% 2. Bed shear stress (Shields parameter) for Trough half cycle
%    alpha VA2013 Eqn. 19
%-----------------------------------------------------------------------
%
[eta_t, dsf_t, mag_theta_t, tau_wRe_t] = ...
    stress_progressive_surface_waves(d50, d90, osmgd,         ...
    Td, depth,                ...
    umag_curr, uhat, uhat_t, ahat, ...
    T_tu, T_t, mag_ut);
%
%-----------------------------------------------------------------------
% Shields parameter in trough cycle VA2013 Eqn 15
% rho0 is required for non-dimensionalizing
%-----------------------------------------------------------------------
%
% CRS - I think the tau_wRe term should be added to a negative stress in
% the trough...changed some signs above to propagate a negative ut here.
% CRS - I think the current roughness and stress should be for the whole
% wave cycle, not divided into crest/trough components
% CRS - Santoss has a coefficient in front of the progressive wave
% streaming term eps_eff, and might use a value of 0.55...or zero to negate
% the term
theta_tx=abs(mag_theta_t)*ut_rx/(mag_ut)+eps_eff*tau_wRe_t*osmgd/rho0;
theta_ty=abs(mag_theta_t)*ut_ry/(mag_ut);
mag_theta_t=sqrt(theta_tx*theta_tx+theta_ty*theta_ty);
fprintf(1,'|theta_t|,theta_tx, theta_ty, theta_tprog: %f %f %f %f\n',...
    mag_theta_t, theta_tx, theta_ty, eps_eff*tau_wRe_t*osmgd/rho0 )
%%
%-----------------------------------------------------------------------
% 3. Compute sediment load entrained during each crest half cycle
%-----------------------------------------------------------------------
%
%-----------------------------------------------------------------------
%      Crest half cycle
%-----------------------------------------------------------------------
%
wavecycle=1.0;
[ om_cc, om_ct ]= sandload_vandera(wavecycle,...
    Hs, Td,  depth, RR,                 ...
    d50, rhos, c_w,                     ...
    eta_c, dsf_c,                           ...
    T_c, T_cu, uhat_c, mag_theta_c);
%
%-----------------------------------------------------------------------
%       Trough half cycle
%-----------------------------------------------------------------------
%
wavecycle=-1.0;
[om_tt, om_tc] = sandload_vandera(wavecycle,...
    Hs, Td,  depth, RR,                 ...
    d50, rhos, c_w,                     ...
    eta_t, dsf_t,                           ...
    T_t, T_tu, uhat_t, mag_theta_t);

%
%-----------------------------------------------------------------------
% VA2013  Use the velocity-load equation 1.
% Non-dimensional net transport rate
%-----------------------------------------------------------------------
%
smgd_3=sqrt((rhos/rho0-1.0)*g*d50.^3.0);
%
cff1=0.5*T_c/(T_cu);
cff2=sqrt(mag_theta_c)*T_c*(om_cc+cff1*om_tc);
%
cff3=theta_cx/mag_theta_c;
bedld_cx=cff2*cff3;
%
cff3=theta_cy/mag_theta_c;
bedld_cy=cff2*cff3;

cff1=0.5*T_t/(T_tu);
cff2=sqrt(mag_theta_t)*T_t*(om_tt+cff1*om_ct);

cff3=theta_tx/mag_theta_t;
bedld_tx=cff2*cff3;
%
cff3=theta_ty/mag_theta_t;
bedld_ty=cff2*cff3;
%
% The units of these are m2 sec-1
% bed_frac, rhos multiplied
%
bed_frac=1.0;
bedld_x=bed_frac*rhos*smgd_3*(bedld_cx+bedld_tx)/Td;
bedld_y=bed_frac*rhos*smgd_3*(bedld_cy+bedld_ty)/Td;


%% sandload_vanderaa
function [om_ii, om_iy]=sandload_vandera(wavecycle,...
    Hs, Td,  depth, RR,                 ...
    d50, rhos, c_w,                     ...
    eta, dsf,                           ...
    T_i, T_iu, uhat_i, mag_theta_i)
%
% VA2013 Text under equation 37
%
m=11.0; n=1.2; alpha=8.2;
xi=1.7; % Based on Santoss_core.m
%
%
% Find settling velocity based on Soulsby (1997).
% VA2013 Use 0.8*d50 for settling velocity (text under equation 28).
%
w_s=w_s_calc(0.8*d50, rhos);
%
% VA2013 Equation 29, for crest cycle
%
ws_eta=w_sc_calc(Hs, Td, depth, RR, w_s, eta);
ws_dsf=w_sc_calc(Hs, Td, depth, RR, w_s, dsf);
if(wavecycle==1.0)
    w_sc_eta=max(w_s+ws_eta,0.0);
    w_sc_dsf=max(w_s+ws_dsf,0.0);
end
%
% VA2013 Equation 30, for trough cycle
% Upward orbital velocity decreases ws.
% VA2013 allows this to go to zero. I think that is ok, but Taran put a
% limit on it at 0.36*ws, citing
% Note that VA2013 uses 0.8*D50 for suspended sediment...we are not doing
% that
%
if(wavecycle==-1.0)
    w_sc_eta=max(w_s-ws_eta,0.36*w_s);
    w_sc_dsf=max(w_s-ws_dsf,0.36*w_s);
end
%
% VA2013 Equation 33, Phase lag parameter
%
cff=(1.0-wavecycle*xi*uhat_i)/c_w
cff1_eta=(1.0/(2.0*(T_i-T_iu)*w_sc_eta));
cff1_dsf=(1.0/(2.0*(T_i-T_iu)*w_sc_dsf));
%        cff1_eta=(1.0/(2.0*(T_iu)*w_sc_eta))
%        cff1_dsf=(1.0/(2.0*(T_iu)*w_sc_dsf))
%
% For ripple regime
%
% TODO: this can be made a little more robust
% if(eta>0.0)
%     P=alpha*eta*cff*cff1_eta;
% elseif(eta==0.0)
%     %
%     % For sheet flow regime
%     %
%     P=alpha*dsf*cff*cff1_dsf;
% end
% CRS like this:
P=alpha*dsf*cff*cff1_dsf;
if(eta>d50)
    P=alpha*eta*cff*cff1_eta;
end

% TODO CRS - this is not the same eps_eff that I put into the main routine
eps_eff=(dsf/d50)^0.25; 
theta_ieff=eps_eff*mag_theta_i;
%
% Find critical Shields parameters based on Soulsby (1997).
%
theta_cr=theta_cr_calc(d50, rhos);
%
% Sand load entrained in the flow during each half-cycle
%
om_i=max( (theta_ieff-theta_cr), 0.0);
om_i=m*(theta_ieff-theta_cr)^n;
%
% VA2013 Equation 23-26, Sandload entrained during half cycle
%
if(P<=1.0)
    om_ii=om_i;
    om_iy=0.0;
else
    om_ii=om_i/P;
    cff=1.0/P;
    om_iy=om_i*(1.0-cff);
end
%
return
end % function sandload_vandera

%% w_sc_calc
%
function worb = w_sc_calc(Hs, Td, depth, RR, w_s, zws)
%
% Second order Stokes theory to get vertical velocity of water particle
% at a given elevation based on santoss_core.m
%
worb1=pi*Hs*zws/(Td*depth);
worb2=worb1*2.0*(RR+RR-1.0);
%
%  Using the SANTOSS model formulation
%
cff=1.0/8.0;
worb=cff*worb1*sqrt(64.0-(-worb1+...
    sqrt(worb1^2.0+32.0*...
    worb2^2.0))^2.0/(worb2^2.0))+...
    worb2*sin(2.0*acos(cff*(-worb1+...
    sqrt(worb1^2.0+32.0*worb2^2.0))/worb2));

return
end % function w_sc_calc

%% stress_progressive_surface_waves
function [eta, dsf, mag_theta_i, tau_wRe] = ...
    stress_progressive_surface_waves(d50, d90, osmgd,...
    Td, depth,...
    umag_curr, uhat, uhat_i, ahat,...
    T_iu, T_i, mag_ui )
%
% Input the crest or trough half cycle velocity
% d50 -- grain size in meters
% Different for crest and trough half cycles
%
tol = 0.001;
total_iters=10;
rho0 = 1025.0;
%
% Iterative solution to obtain current and wave related bed roughness
% VA2013 Apendix A, Shields parameter (Stress) depends on bed roughness
% Bed roughness computed from converged Shields parameter
%
% maximum mobility number at crest and trough
% For irregular waves, use Rayleigh distributed u(1/10) value
% VA, text under equation Appendix B.4
%
% smgd=(rhos/rho0-1.0)*g*d50
% osmgd=1.0/smgd
psi=(1.27*uhat)^2*osmgd;
%
% Use Appendix B eqn B.1 and B.2 to get ripple height and length
%
[eta, lambda] = ripple_dim(psi, d50);
% These agree:
[eta_CRS, lambda_CRS]=od_ripple(d50, psi);
checkvals('eta',eta,eta_CRS);
checkvals('lambda',lambda,lambda_CRS);
%
eta=eta*ahat;
lambda=lambda*ahat;

% VA2013 Eqn. 19:
alpha=umag_curr/(umag_curr+uhat);
%
% Initiliaze with theta_timeavg=0 and theta_hat_i=theta_timeavg
%
theta_timeavg=0.0;
theta_timeavg_old=0.0;
theta_hat_i=theta_timeavg;
for iter=1:total_iters
    %
    % Sheet flow thickness VA2013 Appendix C C.1
    %
    dsf=dsf_calc(d50, theta_hat_i); %this dsf is in m
%     % These agree
    dsf_CRS = dsf_func(d50, theta_hat_i);
    checkvals('dsf',dsf, dsf_CRS)
    %
    % Calculate current-related bed roughness from VA2013 Appendix A.1
    %
    ksd=ksd_calc(d50, d90, mu_calc(d50), theta_timeavg, eta, lambda);
    %
    % Calculate full-cycle current friction factor from VA2013 Eqn. 20
    %
    fd=fd_calc(dsf, ksd);
    fd_CRS = fd_func(dsf, ksd);
    checkvals('fd',fd, fd_CRS);
    %
    % Calculate wave related bed roughness from VA2013 A.5
    %
    ksw=ksw_calc(d50, mu_calc(d50), theta_timeavg, eta, lambda);
    ksw_CRS = ksw_func(d50, eta, lambda, theta_timeavg);
    checkvals('ksw',ksw, ksw_CRS);
    %
    % Calculate full-cycle wave friction factor VA2013 Appendix Eqn. A.4
    %
    fw=fw_calc(ahat, ksw);
    fw_CRS = fw_func(ahat, ksw);
    checkvals('fw',fw, fw_CRS);
    %
    %
    % Calculate Time-averaged absolute Shields stress VA2013 Appendix Eq. A.3
    %
    theta_timeavg=osmgd*(0.5*fd*umag_curr^2.0+...
        0.25*fw*uhat^2.0);
    %
    % Wave friction factor for wave and crest half cycle VA2013 Eqn. 21
    %
    fw_i=fwi_calc(T_iu, T_i, ahat, ksw);
    %
    % Wave current friction factor (Madsen and Grant) VA2013 Eqn. 18
    % Different for crest and trough
    %
    fwd_i=alpha*fd+(1.0-alpha)*fw_i;
    %
    % Update theta_hat_i based on crest/trough amplitude uhat Eqn. C.2
    %
    theta_hat_i=0.5*fwd_i*uhat_i^2*osmgd;
    %
    fprintf(1,'%d %f %f %f %f %f\n',iter,dsf,fd,fw,theta_timeavg,theta_hat_i)
    if(abs(theta_timeavg-theta_timeavg_old) < tol)
        break
    end
    if((abs(theta_timeavg-theta_timeavg_old) >= tol)&&(iter==total_iters))
        fprintf(1,'Warning...stress calcs did not converge.\n')
    end
    theta_timeavg_old=theta_timeavg;
end

%
% Recompute sheet flow thickness to use converged values of theta_timeavg
% and theta_hat_i
%
% Sheet flow thickness VA2013 Appendix C C.1
% Update from initial value
%
dsf=dsf_calc(d50, theta_hat_i); %this dsf is in m
%
% Calculate current-related bed roughness from VA2013 Appendix A.1
%
ksd=ksd_calc(d50, d90, mu_calc(d50), theta_timeavg, eta, lambda);
%
% Calculate full-cycle current friction factor from VA2013 Eqn. 20
%
fd=fd_calc(dsf, ksd);
%
% Calculate wave related bed roughness from VA2013 A.5
%
ksw=ksw_calc(d50, mu_calc(d50), theta_timeavg, eta, lambda);
%
% Calculate full-cycle wave friction factor VA2013 Appendix Eqn. A.4
%
fw=fw_calc(ahat, ksw);
%
% Calculate wave Reynolds stress from full cycle wave and friction factor
% that were formed from the iterative cycle, VA2013, Eqn.22
%
fwd=alpha*fd+(1.0-alpha)*fw;
%
k=kh_calc(Td,depth)/depth;     % Wave number
c_w=2*pi/(k*Td);               % Wave speed
alpha_w=0.424;
%
tau_wRe=rho0*fwd*alpha_w*uhat^3.0/(2.0*c_w);
fprintf(1,'------------------------------\n')
fprintf(1,'tau_wRe %f\n',tau_wRe)
%
% Wave friction factor for wave and crest half cycle VA2013 Eqn. 21
%
fw_i=fwi_calc(T_iu, T_i, ahat, ksw);
%
% Wave current friction factor (Madsen and Grant) VA2013 Eqn. 18
% Different for crest and trough
%
fwd_i=alpha*fd+(1.0-alpha)*fw_i;
%
% VA2013-Magnitude of Shields parameter Eqn. 17
%
mag_theta_i=0.5*fwd_i*mag_ui^2*osmgd;
end % function stress_progressive_surface_waves

%% ripple_dim
function [eta, lambda]=ripple_dim(psi, d50)
%
% Calculate ripple dimensions of O'Donoghue et al. 2006
% based on VA2013 Appendix B
% returns eta-Ripple length and lambda-Ripple length
d50_mm=0.001*d50;
if(d50_mm<0.22)
    m_eta=0.55;
    m_lambda=0.73;
elseif(d50_mm>0.22 && d50_mm<0.30)
    m_eta=0.55+(0.45*(d50_mm-0.22)/(0.30-0.22));
    m_lambda=0.73+(0.27*(d50_mm-0.22)/(0.30-0.22));
else
    m_eta=1.0;
    m_lambda=1.0;
end
%
% Smooth transition between ripple regime and bed sheet flow regime
%
if(psi<=190.0)
    n_eta=1.0;
elseif(psi>190.0 && psi<240.0)
    n_eta=0.5*(1.0+cos(pi*(psi-190.0)/(50.0)));
elseif(psi >= 240.0)
    n_eta=0.0;
end
n_lambda=n_eta;
%
eta=max(0.0,m_eta*n_eta*(0.275-0.022*psi^0.42));
lambda=max(0.0,m_lambda*n_lambda*(1.97-0.44*psi^0.21));
%
return
end % function ripple_dim

%% theta_cr_calc        %
function theta_cr = theta_cr_calc(d50, rhos)
%
% Critical Shields parameter from Soulsby (1997).
%
rho0=1025.0;
nu = 1.36E-6;
g = 9.81;

s=rhos/rho0;
dstar=(g*(s-1)/(nu*nu))^(1.0/3.0)*d50;
cff1=0.30/(1.0+1.2*dstar);
cff2=0.055*(1.0-exp(-0.020*dstar));
theta_cr=cff1+cff2;
%
return
end % function theta_cr_calc


%% w_s_calc
function ws = w_s_calc(d50, rhos)
%
% Critical Shields parameter from Soulsby (1997).
% Dynamics of Marine Sands
rho0=1025.0;
nu = 1.36E-6;
g = 9.81;

s=rhos/rho0;
dstar=(g*(s-1)/(nu*nu))^(1.0/3.0)*d50;
cff=nu/d50;
cff1=10.36;
ws=cff*(sqrt(cff1*cff1+1.049*dstar^3.0)-cff1);
%
return
end % function w_s_calc


%% mu_calc
function mu = mu_calc(d50)
%
% Calculate bed roughness factor based on grain size
% VA2013 Appendix A., required for current related bed roughness
% and wave related bed roughness.
d50_mm=d50*0.001;
%
if(d50_mm<=0.15)
    mu=6.0;
elseif(d50_mm>0.15 && d50_mm<0.20)
    mu=6.0-5.0*((d50_mm-0.15)/(0.2-0.15));
elseif(d50_mm>0.20)
    mu=1.0;
end
return
end % function mu_calc


%% ksd_calc
function ksd = ksd_calc(d50, d90, mu, theta_timeavg, eta, lambda)
%
% Calculate current-related bed roughness from VA2013 Appendix A.1.
%

%
eta=max(eta,d50);
lambda=max(lambda,d50);
ripple_fac=0.4*eta^2.0/lambda;
ksd=max( 3.0*d90,...
    d50*(mu+6.0*(theta_timeavg-1.0)) )+...
    ripple_fac;
%
return
end % function ksd_calc

%% ksw_calc
function ksw = ksw_calc(d50, mu, theta_timeavg, eta, lambda)
%
% Calculate wave related bed roughness from VA2013 Eqn. A.5.
%
eta=max(eta,d50);
lambda=max(lambda,d50);
ripple_fac=0.4*eta^2.0/lambda;
ksw = max( d50,...
    d50*(mu+6.0*(theta_timeavg-1.0)) )...
    +ripple_fac;
%
return
end % function ksw_calc

%% fw_calc
function fw = fw_calc(ahat, ksw)
%
% Calculate full-cycle wave friction factor from VA2013 Eqn. A.4.
%
ratio=ahat/ksw;
if(ratio>1.587)
    fw = 0.00251*exp(5.21*(ratio)^(-0.19));
else
    fw = 0.3;
end
%
return
end % function fw_calc

%% fd_calc
function fd = fd_calc(dsf, ksd)
%
% Calculate current related friction factor VA2013 Eqn. 20
% Assuming logarithmic velocity profile.
%
von_k = 0.41;
fd = 2.0*(von_k/log(30.0*dsf/ksd))^2.0;
%
return
end % function fd_calc


%% fwi_calc
function fwi = fwi_calc(T_iu, T_i, ahat, ksw)
%
% Wave friction factor for wave and crest half cycle VA2013 Eqn. 21.
%
c1=2.6;
ratio=ahat/ksw;
if(ratio>1.587)
    cff=(2.0*T_iu/T_i)^c1;
    fwi = 0.00251*exp(5.21*(cff*ratio)^(-0.19));
else
    fwi = 0.3;
end
%
return
end % function fwi_calc


%% dsf_calc
function dsf = dsf_calc(d50, theta_i)
%
% Sheet flow thickness VA2013 Appendix C.1.
%
d50_mm=d50*0.001;
if(d50_mm<=0.15)
    cff=25.0*theta_i;
elseif(d50_mm>0.15 & d50_mm<0.20)
    cff=25.0-(12.0*(d50_mm-0.15)/0.05);
elseif(d50_mm >= 0.20)
    cff=13.0*theta_i;
end
dsf = max(d50*cff,d50);
%
return
end % function dsf_calc
%
%  end of functions for step 2 for Shear stress formulation
%

%% skewness_params
function [r, phi] = skewness_params(H_s, T, depth)
%
% Ruessink et al. provides equations for calculating skewness parameters
% Uses Malarkey and Davies equations to get "bb" and "r"
% Given input of H_s, T and depth
% r     - skewness/asymmetry parameter r=2b/(1+b^2)            [value]
% phi   - skewness/asymmetry parameter                         [value]
% Su     - umax/(umax-umin)                                    [value]
% Au   - amax/(amax-amin)                                      [value]
% alpha - tmax/pi                                              [value]

p1=0.0;
p2=0.857;
p3=-0.471;
p4=0.297;
p5=0.815;
p6=0.672;
dtr = pi/180.
%
% Ruessink et al., 2012, Coastal Engineering 65:56-63.
%
% k is the local wave number computed with linear wave theory.
%
k=kh_calc(T,depth)/depth;
%
% H_s=sqrt(2.0)*H_rms
a_w=0.5*H_s;
Ur=0.75*a_w*k/((k*depth)^3.0);
%
% Ruessink et al., 2012 Equation 9.
%
cff=exp( (p3-log10(Ur)) /p4);
B=p1+((p2-p1)/(1.0+cff));
psi=-90.0*dtr*(1.0-tanh(p5/Ur^p6));
%
% Markaley and Davies, equation provides bb which is "b" in paper
% Check from where CRS found these equations
%
bb=sqrt(2.0)*B/(sqrt(2.0*B^2.0+9.0));
r=2.0*bb/(bb^2.0+1.0);
%
% Ruessink et al., 2012 under Equation 12.
%
phi=-psi-0.5*pi;
%
% Where are these asymmetry Su, Au utilized
% recreate the asymetry
%
Su=B*cos(psi);
Au=B*sin(psi);

return
end % function skewness_params

%% abreu_points
function [DTc, DTt, DTcu, DTtu, umax, umin, RR] = ...
    abreu_points(r, phi, Uw, T )
%
%  Calculate umax, umin, and phases of asymmetrical wave orbital velocity
%
%  Use the asymmetry parameters from Ruessink et al, 2012
%  to get the umax, umin and phases of asymettrical wave
%  orbital velocity to be used by Van Der A.
%  T_c is duration of crest
%  T_cu Duration of accerating flow within crest half cycle
%
omega=2.0*pi/T;
%
phi_new=-phi;

% Malarkey and Davies (Under equation 16b)
P=sqrt(1.0-r*r);
%
% Malarkey and Davies (Under equation 16b)
%
b=r/(1.0+P);
%
% Appendix E of Malarkey and Davies
%
c=b*sin(phi_new);
%
% Appendix E of Malarkey & Davies
% Phase of umax (crest) and umin (trough) (in radians, from 0 to 2*pi)
c = b*sin(phi);
ratio_c = ((4.*c*(b*b-c*c)+(1.-b*b)*(1.+b*b-2.*c*c))/((1.+b*b).^2-4.*c*c));
ratio_c = min(ratio_c,1.0);
ratio_c = max(ratio_c,-1.0);
tmc = asin(ratio_c);
if(tmc<0.)
   tmc = tmc+2*pi;
end
ratio_t = ((4.*c*(b*b-c*c)-(1.-b*b)*(1.+b*b-2.*c*c))/((1.+b*b).^2-4.*c*c));
ratio_t = min(ratio_t,1.0);
ratio_t = max(ratio_t,-1.0);
tmt = asin(ratio_t);
if(tmt<0.)
   tmt = tmt+2.*pi;
end
%
% Non dimensional umax and umin, under E5 in Malarkey and Davies
%
umax=1.0+c;
umin=umax-2.0;
%
%       Dimensionalize
%
umax=umax*Uw;
umin=umin*Uw;
%
% phase of zero upcrossing and downcrossing (radians)
%
tzu=asin(b*sin(phi_new));
tzd=2.0*acos(c)+tzu;
%
% MD, equation 17
%
RR=0.5*(1.0+b*sin(phi_new));
%
% MD, under equation 18
%
if(r<=0.5)
    F0=1.0-0.27*(2.0*r)^(2.1);
else
    F0=0.59+0.14*(2.0*r)^(-6.2);
end
%
% MD, Equation 15a,b
%
if(r >= 0.0 && r<0.5)
    betar_0=0.5*(1.0+r);
elseif(r>0.5 && r<1.0)
    cff1=4.0*r*(1.0+r);
    cff2=cff1+1.0;
    betar_0=cff1/cff2;
end
%
% MD, Equation 18
%
cff=sin((0.5*pi-abs(phi_new))*F0)/sin(0.5*pi*F0);
beta=0.5+(betar_0-0.5)*cff;
%
% MD, Table 1, get asymmetry parameterization
% using GSSO (10a,b)
%
cff=sqrt(2.0*(1.0+b*b)^3.0);
Sk=3.0*sin(phi_new)/cff;
Ak=-3.0*cos(phi_new)/cff;
%
% These are the dimensional fractions of wave periods needed by Van der A eqn.
% TSK - Check source of these equations
%
w=1.0/omega;
DTc=(tzd-tzu)*w;
DTt=T-DTc;
DTcu=(tmc-tzu)*w;
DTtu=(tmt-tzd)*w;
%
T_tu=tzd*w;
T_cu=tzu*w;
T_c=tmc*w;
T_t=tmt*w;
%
return
end % function abreu_points

%% kh_calc
function kh = kh_calc(Td,depth)
%
%  Calculate wave number from Wave period and depth
%
% RL Soulsby (2006) "Simplified calculation of wave orbital velocities"
% HR Wallingford Report TR 155, February 2006
%
g = 9.81;
omega=2.0*pi/Td;
x=omega^2.0*depth/g;
%
if(x<1.0)
    y=sqrt(x);
else
    y=x;
end
%
% Iterating three times to solve eqn.7 of Soulsby 1997 by using
% eqns. (12a-14)
%
t=tanh(y);
cff=(y*t-x)/(t+y*(1.0-t*t));
y=y-cff;
%
t=tanh(y);
cff=(y*t-x)/(t+y*(1.0-t*t));
y=y-cff;
%
t=tanh(y);
cff=(y*t-x)/(t+y*(1.0-t*t));
y=y-cff;
kh=y;
%
return
end % function kh_calc

%% checkvals
function checkvals(name, val1, val2, tol)
if (exist('tol')~=1),tol = 1e-4; end
if( abs(val1-val2)>tol)
    fprintf('WARNING: Diff. values for %s: %f %f\n',name,val1,val2);
end
return
end
