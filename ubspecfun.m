function [ub,Tbav]=ubspecfun(hs,tp,h)
% UBSPECFUN - Caclulate ub and Tbav from Hs and Tp at surface using JONSWAP
% [ub,Tbav] = ubspecfun(hs,tp,h)
%
% Input:
%   hs - Significant wave height (m)
%   tp - Peak period (s)
%   h  - Water depth (m)
%
% Written by Pat Wiberg, UVa
% Minor changes by csherwood@usgs.gov
% Not sure if 
g=9.81;
f=30:10:400; f=f./1000; sf=2.*pi.*f; nb=length(f); Tbin=1./f;

ub=zeros(length(hs),1); Hb=ub; Tbav=ub;
% for i=1:length(tp),
%     khd(i)=rtnewt((2*pi./tp(i)).^2*h./g);
% end;
khd=qkhf( 2*pi/tp, h );
khd=khd(:);
ubTd=(pi./tp).*hs./sinh(khd);

sjmat=[];
for i=1:length(hs),
    fp=1./tp(i); sfp=2.*pi*fp;
    gam=1.25; sig=0.08; alp=0.0081;
    if fp~=Inf,
        eterm=-((sf-sfp).^2)/(2.*sig.^2.*sfp.^2);
        ee=exp(eterm);
        t2=gam.^ee;
        t1=-1.25.*(sfp./sf).^4;
        sjmat(i,:)=alp*g.^2./sf.^5.*exp(t1).*t2;
    else,
        sjmat(i,:)=zeros(1,nb);
    end
end;

cfmat=[];
%for i=1:10,
for i=1:length(tp),
    ubT=zeros(1,nb);
    for k=1:nb,
        kh(k)=qkhf( 2*pi/Tbin(k), h );
        %kh(k)=rtnewt((2*pi./Tbin(k)).^2.*h./g);
    end;
    iter=0; cf=2;
    while cf<0.99|cf>1.01;
        ubT=(pi./Tbin.*2.*sqrt(sjmat(i,:).*0.01)./sinh(kh)).^2;
        HbT=4.*sjmat(i,:).*.01./((cosh(kh)).^2);
        waT=4.*sjmat(i,:).*0.01;
        ub(i)=sqrt(2).*sqrt(sum(ubT));
        Hb(i)=2.*sqrt(sum(HbT));
        Hss(i)=2.*sqrt(sum(waT));
        fbav(i)=sum(ubT./(Tbin.^2))./sum(ubT);
        Tbav(i)=1./sqrt(fbav(i));
        cf=(hs(i)./Hss(i));
        iter=iter+1;
        if iter>10, break; end;
        sjmat(i,:)=sjmat(i,:).*cf;
    end
    cfmat=[cfmat; cf iter];
end;

