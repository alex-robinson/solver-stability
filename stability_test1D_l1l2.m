function [fac u0 Lm] = stability_test1D_l1l2(...
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

u0 = rhog * H0 / beta_sl * alpha;

N = round(L/delx);
delx = L/N;

if (Lm > 0.2 * L);
    error('membrane stress too long for domain');
end

h = epsh * randn(N,1);
h = [h; h(1)];

Hinit = h + H0;
H = Hinit;

for i = 1:nt;

 rhs = rhog * (H(1:end-1)+H(2:end))/2 .* (diff(H)/delx - alpha);
 %F2 = (H(1:end-1)+H(2:end))/2/(3*eta);
 B = beta_sl*ones(size(rhs));

 d0 = -B - 4*eta*H(1:end-1)/delx^2 - 4*eta*H(2:end)/delx^2;
 dr = 4 * eta * H(2:end) / delx^2;
 dl = 4 * eta * H(1:end-1) / delx^2;
 A = spdiags([dr d0 dl],[-1 0 1],N,N)'; 
 A(1,N) = dl(1);
 A(N,1) = dr(end);

 u = A\rhs;
 u_ext = [u(end); u; u(1)];
 H_ext = [H(end); H];
 H_ext2 = [H(end); H; H(1)];

 Hmid = (H(1:end-1)+H(2:end))/2;
 uav = u + ...
     Hmid.^2/3 .* (4/delx^2 * (-2*u_ext(2:end-1)+u_ext(1:end-2)+u_ext(3:end)) - rhog/eta*(diff(H)/delx - alpha)) + ...
     2 * Hmid .* (diff(H)/delx - alpha) .* (-u_ext(1:end-2)+u_ext(3:end))/2/delx;
 uav_ext = [uav(end); uav; uav(1)];
 flux = (uav_ext .* H_ext);
 
 H = H + del_t/delx * (flux(1:end-1)-flux(2:end));
%     del_t * rhog / (3*eta*delx^2) * H0^3 * (-2*H_ext2(2:end-1).^1 + H_ext2(3:end).^1 + H_ext2(1:end-2).^1);
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
