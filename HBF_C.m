function [F, R_vec, sigma_max] = HBF_C(H,Nr,Nt,N_vec,nbar,L,K,q,rho)
sigma_max = 0;
F = zeros(L*Nr,L*nbar);
I = eye(K*Nt);
Q = I;
R_vec = zeros(L,1);
col_end = 0;
for l = 1:L
    N_l = N_vec(l);
    col_start = col_end + 1;
    col_end = col_start + N_l - 1;
    
    H_l = H((l-1)*Nr+1:(l-1)*Nr+Nr, :);
    T_l = H_l*(Q^(-1))*H_l';
    [U,D,V] = svd(T_l);
    F_l = U(:,1:N_l);
    %F_l*F_l'
    F_l = quant_sub(Nr,L,N_l,F_l,q); % quantize U
    R_vec(l) = log2(det(eye(N_l) + rho*(F_l' * T_l * F_l)));
    
    %sigma_tmp = diag(D); sigma = sigma_tmp(1:N_l); sigma_max(l) = sigma(1);
    %R_vec1(l) = sum(log2(1 + rho*sigma));
    
    G_l = H_l'*F_l*F_l'*H_l;
    E_l = I + rho*Q^(-1)*G_l;
    Q = Q * E_l;
    F((l-1)*Nr+1:(l-1)*Nr+Nr, col_start:col_end) = F_l;
end

if (size(F,2) ~= L*nbar)
    error('***********');
end
end