function [anl] = analtical_hybrid(dx,alpha,mu,H0,beta_sl)
    % Calculate analytical solution - hybrid
    
    rhog = 910*9.8;

    F2 = H0/(3*mu);
    B  = beta_sl;
    D  = 4*mu*H0./dx.^2;
    Q  = B./D;
    r  = 1 + Q./2 - sqrt((1+Q/2).^2-1);
    a0 = 1./(2+Q-2.*r);

    u0   = rhog * H0 * alpha * (1/beta_sl + H0/(3*mu));
    ueff = rhog * H0 * alpha * (1/beta_sl + H0/(  mu));
    
    dt_adv   = dx/ueff;
     
    % Numerically test dt values against stability critierion
    dt_dyn = dx*NaN;
    for i=1:length(dt_dyn);
        dt_dyn(i) = analtical_hybridl1l2_iter(dx(i),alpha,mu,H0,beta_sl,"hybrid");
    end
    
    % Populate anl struct for output:
    anl.dx = dx;
    anl.u0 = u0;
    anl.dt_adv   = dt_adv;
    anl.dt_dyn   = dt_dyn;
    
end
