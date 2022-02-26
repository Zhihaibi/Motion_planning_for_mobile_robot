function [Aeq beq]= getAbeq(n_seg, n_order, waypoints, ts, start_cond, end_cond)
    n_all_poly = n_seg*(n_order+1);
    %#####################################################
    % p,v,a,j constraint in start, 
    Aeq_start = zeros(4, n_all_poly);
    beq_start = zeros(4, 1);
    % STEP 2.1: write expression of Aeq_start and beq_start
    % v,a,j = 0，First segment (t = 0)
    Aeq_start(1,1) = 1; % p
    Aeq_start(2,2) = 1; % v
    Aeq_start(3,3) = 2; % a
    Aeq_start(4,4) = 6; % J
    beq_start(1) = start_cond(1); % p
    
    %#####################################################
    % p,v,a,j constraint in end
    Aeq_end = zeros(4, n_all_poly);
    beq_end = zeros(4, 1);
    % STEP 2.2: write expression of Aeq_end and beq_end
    % v,a,j = 0,Last segment (t = ts)
    Aeq_end(1,end-7:end)  = flip([ts(n_seg)^7, ts(n_seg)^6, ts(n_seg)^5, ts(n_seg)^4, ts(n_seg)^3, ts(n_seg)^2, ts(n_seg), 1]);     % p
    Aeq_end(2,end-7:end)  = flip([7*ts(n_seg)^6, 6*ts(n_seg)^5, 5*ts(n_seg)^4, 4*ts(n_seg)^3, 3*ts(n_seg)^2, 2*ts(n_seg)^1, 1, 0]); % v
    Aeq_end(3,end-7:end)  = flip([7*6*ts(n_seg)^5, 6*5*ts(n_seg)^4, 5*4*ts(n_seg)^3, 4*3*ts(n_seg)^2, 3*2*ts(n_seg)^1, 2, 0 , 0]);  % a
    Aeq_end(4,end-7:end)  = flip([7*6*5*ts(n_seg)^4, 6*5*4*ts(n_seg)^3, 5*4*3*ts(n_seg)^2, 4*3*2*ts(n_seg)^1, 3*2*1, 0, 0, 0]);     % J
    beq_end(1) = end_cond(1); % P
    
    %#####################################################
    % position constrain in all middle waypoints
    Aeq_wp = zeros(n_seg-1, n_all_poly);
    beq_wp = zeros(n_seg-1, 1);
    % STEP 2.3: write expression of Aeq_wp and beq_wp
    % p1 = p2，中间点的位移，第二段曲线0时刻至最后一段曲线0时刻的值
    for i = 1:n_seg-1
        Aeq_wp(i,i*8+1) = 1;
        beq_wp(i) = waypoints(i+1);
    end
    
    %#####################################################
    % position continuity constrain between each 2 segments
    Aeq_con_p = zeros(n_seg-1, n_all_poly);
    beq_con_p = zeros(n_seg-1, 1);
    % STEP 2.4: write expression of Aeq_con_p and beq_con_p
    % p1(t) = p2(0)，两段曲线之间的连接，前一段曲线ts时刻的位移等于后一段曲线0时刻的位移

    %#####################################################
    % velocity continuity constrain between each 2 segments
    Aeq_con_v = zeros(n_seg-1, n_all_poly);
    beq_con_v = zeros(n_seg-1, 1);
    % STEP 2.5: write expression of Aeq_con_v and beq_con_v
    % V1(t) = V2(0)

    %#####################################################
    % acceleration continuity constrain between each 2 segments
    Aeq_con_a = zeros(n_seg-1, n_all_poly);
    beq_con_a = zeros(n_seg-1, 1);
    % STEP 2.6: write expression of Aeq_con_a and beq_con_a
    % a1(t) = a2(0)

    %#####################################################
    % jerk continuity constrain between each 2 segments
    Aeq_con_j = zeros(n_seg-1, n_all_poly);
    beq_con_j = zeros(n_seg-1, 1);
    % STEP 2.7: write expression of Aeq_con_j and beq_con_j
    % J1(t) = J2(0)
 
    %#####################################################
    % combine all components to form Aeq and beq   
    Aeq_con = [Aeq_con_p; Aeq_con_v; Aeq_con_a; Aeq_con_j];
    beq_con = [beq_con_p; beq_con_v; beq_con_a; beq_con_j];
    
    for k = 0:1:n_seg-2 % here k is the index of segments
        Aeq_con(1+4*k:4+4*k,1+8*k:8+8*k) = getCoeff(ts(k+1));
        Aeq_con(1+4*k:4+4*k,1+8*(k+1):8+8*(k+1)) = -getCoeff(0);            
    end
    
    
    Aeq = [Aeq_start; Aeq_end; Aeq_wp; Aeq_con];
    beq = [beq_start; beq_end; beq_wp; beq_con];
end