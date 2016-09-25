function gamma = coherence(lam,delta_lam,tau)
 gamma = exp( -(tau*2*pi*delta_lam/lam^2)^2/2 );
end