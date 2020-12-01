function [R, comp, N_vec] = ARFA_ES(H,Nr,Nt,N,L,K,q,rho)

R = 0;
N_vals = 1:L*N;
N_all_comb = permn(N_vals, L);
n_comb = size(N_all_comb,1);
N_vec_best = N*ones(L,1);
for n = 1:n_comb
    N_vec = N_all_comb(n,:);
    if (sum(N_vec) == L*N) && (sum(N_vec(N_vec > Nr)) == 0)
        [Ftmp, Rvec] = C_HBF(H,Nr,Nt,N_vec,N,L,K,q,rho);
        Rtmp = log2(det(eye(K*Nt) + rho*(H'*Ftmp * Ftmp'*H)));
        if Rtmp > R
            F = Ftmp; R = Rtmp; N_vec_best = N_vec; Rvec_best = Rvec;
        end
    end
end
F = C_HBF(H,Nr,Nt,N_vec_best,N,L,K,q,rho);
R = log2(det(eye(K*Nt) + rho*(H'*F * F'*H)));
comp = n_comb;

end % eof