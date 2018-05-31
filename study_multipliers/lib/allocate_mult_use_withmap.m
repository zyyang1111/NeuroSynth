function [mult_queue, mult_out_ind, mult_out_num] = allocate_mult_use_withmap(multipliers, mult_queue, mult_out_ind, mult_out_num, mult_use_array, out_ind_base, Tm, w_map)

for m = 1 : length(mult_use_array)
    mult_use = mult_use_array(m).set;  
    
    % after this loop, multipliers saves the unused multipliers
    for i = 1 : size(mult_use, 1)
        wid = mult_use(i,1); 
        num_wid = mult_use(i,2); 
        
        % multid is the multiplier id to perform the wid multiplication
        multid = w_map(w_map(:,1) == wid, 2); 
        mult_queue(multid) = mult_queue(multid) + num_wid; 
        if mult_out_ind(multid, mult_out_ind(multid, 1)+1) == m + out_ind_base
            mult_out_num(multid, mult_out_num(multid,1)+1) = mult_out_num(multid, mult_out_num(multid,1)+1) + num_wid;
        else
            mult_out_ind(multid, 1) = mult_out_ind(multid, 1) + 1;
            mult_out_num(multid, 1) = mult_out_num(multid, 1) + 1;
            mult_out_ind(multid, mult_out_ind(multid, 1)+1) = m + out_ind_base;
            mult_out_num(multid, mult_out_num(multid, 1)+1) = num_wid;
        end
    end
end
        