function [R, EE, comp, power] = ARFA_SC(H,Nr,Nt,nbar,L,K,q,rho,beta,N)
% SC-HBF scheme
comp = 1;

N_vec_conv = N*ones(L,1);
[Ftmp, eig_tab] = HBF_SC(H,Nr,Nt,N_vec_conv,N,L,K,q);
all_eig = eig_tab(:);

all_eig_sort = sort(all_eig, 'descend');
all_eig_max = all_eig_sort(1:L*nbar);
eig_min = all_eig_max(L*nbar);

for l = 1:L
    eig_l = eig_tab(:,l);
    N_vec(l) = sum(eig_l >= eig_min);
end

% sum(N_vec)
if (sum(N_vec) ~= L*nbar)
    error('***********');
end
[F, R_vec] = HBF_SC(H,Nr,Nt,N_vec,nbar,L,K,q);
R = log2(det(eye(K*Nt) + rho*(H'*F * F'*H)));

% compute EE
n_AP = L-sum(N_vec == 0);
power = get_power(L,N,Nr,n_AP,nbar,'ARFA');
EE = R/power;
end %eof
