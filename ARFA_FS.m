function [R, EE, comp, power] = ARFA_FS(H,Nr,Nt,nbar,L,K,q,rho,N)

N_UB = N;
N_vec = nbar*ones(L,1);
N_vec_best = N_vec;

[Ftmp, R_vec] = HBF_C(H,Nr,Nt,N_vec_best,nbar,L,K,q,rho);
R = log2(det(eye(K*Nt) + rho*(H'*Ftmp * Ftmp'*H)));
[~, ll_order] = sort(R_vec, 'descend');

comp = 1;
while N_vec(ll_order(1)) < N_UB
    tt = 1;
    while (tt < round(L/2))
        kk = L;
        while (N_vec(ll_order(kk)) == 0) && (kk > tt)
            kk = kk - 1;
        end
        if (tt < kk) && (N_vec(ll_order(tt)) < N_UB)
            N_vec(ll_order(tt)) = N_vec(ll_order(tt)) + 1;
            N_vec(ll_order(kk)) = N_vec(ll_order(kk)) - 1;
            if sum(N_vec) == L*nbar
                comp = comp + 1;
                Ftmp = HBF_C(H,Nr,Nt,N_vec,nbar,L,K,q,rho);
                Rtmp = log2(det(eye(K*Nt) + rho*(H'*Ftmp * Ftmp'*H)));
                if Rtmp > R
                    N_vec_best = N_vec; R = Rtmp;
                end
            end
        else
            break
        end
        tt = tt + 1;
    end
end
if (sum(N_vec_best) ~= L*nbar) || (sum(N_vec_best > N_UB) > 0)
    error('***********');
end
F = HBF_C(H,Nr,Nt,N_vec_best,nbar,L,K,q,rho);
R = log2(det(eye(K*Nt) + rho*(H'*F * F'*H)));

% compute EE
n_AP = L-sum(N_vec_best == 0);
power = get_power(L,N,Nr,n_AP,nbar,'ARFA');
EE = R/power;
end % eof