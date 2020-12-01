function [R, EE, comp, power] = ARFA_LS(H,Nr,Nt,nbar,L,K,q,rho,beta,N)
% PL-based SC-ARFA scheme
comp = 1;

N_vec = zeros(L,1);
sum_beta = sum(beta);
for ll = 1:L
    N_vec(ll) = max(round(L*nbar * beta(ll)/sum_beta),1);
    if N_vec(ll) > N
        N_vec(ll) = N;
    end
end
[~, l_order] = sort(beta,'descend');

ii = 1;
while (sum(N_vec) ~= L*nbar)
    if (sum(N_vec) < L*nbar) && (N_vec(l_order(ii)) < N)
        N_vec(l_order(ii)) = N_vec(l_order(ii))+1;
    elseif (sum(N_vec) > L*nbar) && (N_vec(l_order(L-ii+1)) > 0)
        N_vec(l_order(L-ii+1)) = N_vec(l_order(L-ii+1))-1;
    end
    ii = ii+1;
    if ii == L+1
        ii = 1;
    end
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
