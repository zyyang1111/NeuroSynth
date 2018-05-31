function [mult_queue, mult_out_ind, mult_out_num] = allocate_mult_use_withmap_combinelayer(multipliers, mult_queue, mult_out_ind, mult_out_num, mult_use_array, out_ind_base, Tm, w_map)

max_idx = 0;
for m = 1 : length(mult_use_array)
    mult_use = mult_use_array(m).set; 
    if max_idx < size(mult_use, 1)
        max_idx = size(mult_use, 1);
    end
end

for idx = 1 : max_idx
    for m = 1 : length(mult_use_array)
        out_ind_base = mult_use_array(m).base;
        mult_use = mult_use_array(m).set; 
        if idx <= size(mult_use, 1) 
            wid = mult_use(idx, 1);
            num_wid = mult_use(idx, 2); 
        
            % multid is the multiplier id to perform the wid multiplication
            multid = w_map(w_map(:,1) == wid, 2); 
            mult_queue(multid) = mult_queue(multid) + num_wid; 
            if mult_out_ind(multid, mult_out_ind(multid, 1)+1) == out_ind_base
                mult_out_num(multid, mult_out_num(multid,1)+1) = mult_out_num(multid, mult_out_num(multid,1)+1) + num_wid;
            else
                mult_out_ind(multid, 1) = mult_out_ind(multid, 1) + 1;
                mult_out_num(multid, 1) = mult_out_num(multid, 1) + 1;
                mult_out_ind(multid, mult_out_ind(multid, 1)+1) = out_ind_base;
                mult_out_num(multid, mult_out_num(multid, 1)+1) = num_wid;
            end
        end
    end
end
        