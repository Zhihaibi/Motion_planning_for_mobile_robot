function Ct = getCt(n_seg, n_order)
    % STEP 2.1: finish the expression of Ct
    Ct = zeros((n_order+1)*n_seg, 4*n_seg+4);
    k = 0;
    
    for i = 1:(n_order+1)*n_seg/4
        if i == 1
            Ct(1:4,1:4) = eye(4);
            k = k + 4;
            continue;
        end
        
        if i == (n_order+1)*n_seg/4
            k = k + 1;
            Ct((n_order+1)*n_seg-3:(n_order+1)*n_seg, k:k+3) = eye(4);
            continue;
        end
        
        if mod(i,2) == 0
            k = k + 1;
        end
        Ct((i-1)*4+1,k) = 1;
        Ct((i-1)*4+2:(i-1)*4+4, 4*n_seg+2-3*(n_seg-2)+(floor(i/2)-1)*3 : 4*n_seg+4-3*(n_seg-2)+(floor(i/2)-1)*3) = eye(3);
    end
 
end