% test_ruessink_assym
% reproduces Fig. 5 in Ruessink et al., 2012
Ur = logspace( -3, 2., 100)';
Su = zeros( size( Ur ));
Au = zeros( size( Ur ));
ra = ruessink_asymm( Ur );

figure(1); clf
subplot(211)
semilogx( Ur, ra.Su,'linewidth',2)
ylabel('{\itS_u}')
axis([ .3e-2 20., -1 2])
subplot(212)
semilogx( Ur, ra.Au,'linewidth',2)
axis([ .3e-2 20., -2 1])
ylabel('{\itA_u}')
xlabel('{\itU_r}')
