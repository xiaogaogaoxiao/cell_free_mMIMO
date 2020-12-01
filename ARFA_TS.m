function [R, EE, comp, power] = ARFA_TS(H,Nr,Nt,nbar,L,K,q,rho,n_iter,count_max,N)
% TS-aided C-ARFA scheme
comp = 1; count  = 0;

% N = min(L*nbar-(L-1), N);
N_vec = nbar*ones(1,L);
N_vec_best = N_vec;

% initialization
[Ftmp, R_vec] = HBF_C(H,Nr,Nt,N_vec_best,nbar,L,K,q,rho);
R = log2(det(eye(K*Nt) + rho*(H'*Ftmp * Ftmp'*H)));
[~, ll_order] = sort(R_vec, 'descend');
n_larger = sum(R_vec > mean(R_vec));
n_smaller = L - n_larger;

tb_list = N_vec;
ii = 1;
while ii < n_iter
    
    for tt = 1:n_larger
        N_vec_tt = N_vec_best;
        if N_vec_tt(ll_order(tt)) < N
            N_vec_tt(ll_order(tt)) = N_vec_tt(ll_order(tt)) + 1;
        end
        
        for kk = 1:n_smaller
            N_vec_kk = N_vec_tt;
            if N_vec_kk(ll_order(L-kk+1)) > 0
                N_vec_kk(ll_order(L-kk+1)) = N_vec_kk(ll_order(L-kk+1)) - 1;
            end
            if (sum(N_vec_kk) == L*nbar) && (~ismember(N_vec_kk,tb_list,'rows'))
                [Ftmp, R_vec] = HBF_C(H,Nr,Nt,N_vec_kk,nbar,L,K,q,rho);
                Rtmp = log2(det(eye(K*Nt) + rho*(H'*Ftmp * Ftmp'*H)));
                tb_list = cat(1,tb_list,N_vec_kk);
                comp = comp + 1;
                if Rtmp > R
                    N_vec_best = N_vec_kk; R = Rtmp;
                    [~, ll_order] = sort(R_vec, 'descend');
                    n_larger = sum(R_vec > mean(R_vec));
                    n_smaller = L - n_larger;
                    count = 0;
                else
                    count = count + 1;
                end
            end
        end
    end
    if count == count_max
        ii = n_iter;
    end
    ii = ii+1;
end
if (sum(N_vec_best) ~= L*nbar)
    error('***********');
end
[F, R_vec] = HBF_C(H,Nr,Nt,N_vec_best,nbar,L,K,q,rho);
R = log2(det(eye(K*Nt) + rho*(H'*F * F'*H)));

% compute EE
n_AP = L-sum(N_vec_best == 0);
power = get_power(L,N,Nr,n_AP,nbar,'ARFA');
EE = R/power;

end % eof
