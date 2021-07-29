function [fac u0 Lm] = stability_test1D_diva(...
    delx, ...
    del_t, ...
    nt, ...
    beta_sl, ...
    L, ...
    eta, ...
    H0, ...
    alpha, ...
    epsh)

rhog = 910*9.8;

Lm = sqrt(4*eta*H0/beta_sl + 4/3*H0^2);

F2 = H0/(3*eta);
B = beta_sl/(1+beta_sl*F2);
D = 4*eta*H0/delx^2;
Q = B/D;

u0 = rhog * H0 * (1+beta_sl*F2)/beta_sl*alpha;

N = round(L/delx);
delx = L/N;

if (Lm > 0.2 * L);
    error('membranse stress too long for domain');
end

h = epsh * randn(N,1);
h = [h; h(1)];

Hinit = h + H0;
H = Hinit;

for i = 1:nt;

 rhs = rhog * (H(1:end-1)+H(2:end))/2 .* (diff(H)/delx - alpha);
 F2 = (H(1:end-1)+H(2:end))/2/(3*eta);
 B = beta_sl./(1+beta_sl*F2);

 d0 = -B - 4*eta*H(1:end-1)/delx^2 - 4*eta*H(2:end)/delx^2;
 dr = 4 * eta * H(2:end) / delx^2;
 dl = 4 * eta * H(1:end-1) / delx^2;
 A = spdiags([dr d0 dl],[-1 0 1],N,N)'; 
 A(1,N) = dl(1);
 A(N,1) = dr(end);

 u = A\rhs;
 u_ext = [u(end); u; u(1)];
 H_ext = [H(end); H];

 %h = h + u0 * del_t/delx * (h_ext(1:end-1)-h_ext(2:end)) - H0 * del_t/delx * diff(u_ext);
 flux = u_ext .* H_ext;
 H = H + del_t/delx * (flux(1:end-1)-flux(2:end));
 H(end) = H(1);
 
 if ~isempty(find(H<=0));
     break;
 end
 
end

s = std(H);
if (i<100 | isnan(s));
    fac = 1e5;
else
    fac = std(H)/std(Hinit);
end
