function [Aieq, bieq] = getAbieq(n_seg, n_order, corridor_range, ts, v_max, a_max, j_max)
    n_all_poly = n_seg*(n_order+1);
    %#####################################################
    % STEP 3.2.1 p constraint
    Aieq_p = eye(n_all_poly);
    Aieq_p = [Aieq_p; -Aieq_p];
    bieq_p = [];
    
    for i = 1:n_seg
        bieq_p = [bieq_p; corridor_range(i,2)*ones(n_order+1,1)];
    end
    for i = 1:n_seg
        bieq_p = [bieq_p; -corridor_range(i,1)*ones(n_order+1,1)];
    end

    %#####################################################
    % STEP 3.2.2 v constraint   
    Aieq_v_singal = zeros(n_order+1-1,n_order+1);
    for i = 1:1:n_order+1-1
       Aieq_v_singal(i,i) = -n_order;
       Aieq_v_singal(i,i+1) = n_order;
    end
    
    Aieq_v = [];
    for i = 1:n_seg
        Aieq_v = blkdiag(Aieq_v,Aieq_v_singal);
    end
    
    Aieq_v = [Aieq_v;-Aieq_v];
    
    bieq_v = [];
    for i = 1:n_seg*2
        bieq_v = [bieq_v; v_max*ones(n_order+1-1,1)];
    end


    %#####################################################
    % STEP 3.2.3 a constraint   
    Aieq_a_singal = zeros(n_order+1-2,n_order+1);
    for i = 1:n_order+1-2
        Aieq_a_singal(i,i) = n_order*(n_order-1);
        Aieq_a_singal(i,i+1) = -2*n_order*(n_order-1);
        Aieq_a_singal(i,i+2) = n_order*(n_order-1);
    end
    Aieq_a = [];
    for i = 1:n_seg
        Aieq_a = blkdiag(Aieq_a,Aieq_a_singal);
    end
    Aieq_a = [Aieq_a;-Aieq_a];
    
    bieq_a = [];
    for i = 1:n_seg*2
        bieq_a = [bieq_a; a_max*ones(n_order+1-2,1)];
    end

    
    %#####################################################
    % STEP 3.2.4 J constraint   
    Aieq_j_singal = zeros(n_order+1-3,n_order+1);
    for i = 1:n_order+1-3
        Aieq_j_singal(i,i) = n_order * (n_order-1) * (n_order-2);
        Aieq_j_singal(i,i+1) = -3*n_order * (n_order-1) * (n_order-2);
        Aieq_j_singal(i,i+2) = 3*n_order * (n_order-1) * (n_order-2);
        Aieq_j_singal(i,i+3) = -n_order * (n_order-1) * (n_order-2);
    end
    Aieq_j = [];
    for i = 1:n_seg
        Aieq_j = blkdiag(Aieq_j,Aieq_j_singal);
    end
    Aieq_j = [Aieq_j;-Aieq_j];
    
    bieq_j = [];
    for i = 1:n_seg*2
        bieq_j = [bieq_j; j_max*ones(n_order+1-3,1)];
    end
    
    
    %#####################################################
    % combine all components to form Aieq and bieq   
    Aieq = [Aieq_p; Aieq_v; Aieq_a; Aieq_j];
    bieq = [bieq_p; bieq_v; bieq_a; bieq_j];
%     Aieq = Aieq_p;
%     bieq = bieq_p;
end