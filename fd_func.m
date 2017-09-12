function fd = fd_func( dsf, ksd )
% Calculate current friction factor
% fd = fd_func( dsf, ksd )
vk = 0.41;
fd = 2.*(vk/log(30.*dsf/ksd)).^2.; % Eqn. 20
return