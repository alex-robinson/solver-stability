function [anl] = analtical_diva(dx,alpha,mu,H0,beta_sl)
    % Calculate analytical solution - diva
    
    rhog = 910*9.8;

    F2 = H0/(3*mu);
    B  = beta_sl/(1+beta_sl*F2);
    D  = 4*mu*H0./dx.^2;
    Q  = B./D;
    r  = 1 + Q./2 - sqrt((1+Q/2).^2-1);
    a0 = 1./(2+Q-2.*r);

    u0 = rhog * H0 * alpha * (1+beta_sl*F2)/beta_sl;

    dt_adv   = dx/u0;
    dt_dyn   = 2*mu/(rhog*H0)./a0.*(1+r)./(1-r);
    
    % Populate anl struct for output:
    anl.dx = dx;
    anl.u0 = u0;
    anl.dt_adv   = dt_adv;
    anl.dt_dyn   = dt_dyn;
    
end
