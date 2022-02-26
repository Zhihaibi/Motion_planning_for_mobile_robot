function M = getM(n_seg, n_order, ts)
    M = [];
    for k = 1:n_seg
        M_k = zeros(8,8);
        % STEP 1.1: calculate M_k of the k-th segment 
        Coeff1 = getCoeff(0);
        Coeff2 = getCoeff(ts(k));
        M_k = [Coeff1;Coeff2];
        
        M = blkdiag(M, M_k);
    end
end