function [F, eig_tab] = HBF_SC(H,Nr,Nt,N_vec,nbar,L,K,q)

F = zeros(L*Nr,L*nbar);
I = eye(K);
Q = I;
eig_tab = zeros(Nr,L);
col_end = 0;

for l = 1:L
    N_l = N_vec(l);
    col_start = col_end + 1;
    col_end = col_start + N_l - 1;
    
    H_l = H((l-1)*Nr+1:(l-1)*Nr+Nr, :);
    [U,D,V] = svd(H_l);
    
    eigs = diag(D);
    [eig_sort, col_sort] = sort(eigs,'descend');
    eig_tab(1:N_l,l) = eig_sort(1:N_l);
    U = U(:,col_sort(1:N_l));
    F_l = U(:,1:N_l);
    F_l = quant_sub(Nr,L,N_l,F_l,q);

    F((l-1)*Nr+1:(l-1)*Nr+Nr, col_start:col_end) = F_l;
    %R_vec(l) = log2(det(eye(N_l) + rho*(F_l' * H_l * H_l' * F_l)));
end
if (size(F,2) ~= L*nbar)
    error('***********');
end
end