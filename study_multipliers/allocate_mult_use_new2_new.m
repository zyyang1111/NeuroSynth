function [mult_queue, mult_out_ind, mult_out_num] = allocate_mult_use_new2_new(multipliers, mult_queue, mult_out_ind, mult_out_num, mult_use_array, out_ind_base, Tm)

for m = 1 : length(mult_use_array)
    mult_use = mult_use_array(m).set;  
    
    % after this loop, multipliers saves the unused multipliers
    for i = 1 : size(mult_use, 1)
        wid = mult_use(i,1); 
        num_wid = mult_use(i,2); 
        
        dist = 0; 
        multid = find(multipliers(:,1) == wid-dist);
        multid(multipliers(multid, 3) < dist) = []; 
        done_flag = 0; 
        while ~done_flag
            if ~isempty(multid)
                [~, sortid] = sort(multipliers(multid, 2)); % sort according to the variation
                for j = 1 : length(sortid)
                    if multipliers(multid(sortid(j)), 3) >= num_wid
                        % all the operations can be processed using this
                        % type of multipliers
                        mult_queue(multid(sortid(j))) = mult_queue(multid(sortid(j))) + num_wid; 
                        if mult_out_ind(multid(sortid(j)), mult_out_ind(multid(sortid(j)),1)+1) == m + out_ind_base*Tm
                            mult_out_num(multid(sortid(j)), mult_out_num(multid(sortid(j)),1)+1) = ...
                                mult_out_num(multid(sortid(j)), mult_out_num(multid(sortid(j)),1)+1) + num_wid;
                        else
                            mult_out_ind(multid(sortid(j)), 1) = mult_out_ind(multid(sortid(j)), 1) + 1;
                            mult_out_num(multid(sortid(j)), 1) = mult_out_num(multid(sortid(j)), 1) + 1;
                            mult_out_ind(multid(sortid(j)), mult_out_ind(multid(sortid(j)), 1)+1) = m + out_ind_base*Tm;
                            mult_out_num(multid(sortid(j)), mult_out_num(multid(sortid(j)), 1)+1) = num_wid;
                        end

                        num_wid = 0;
                        done_flag = 1; 
                        break;
                        % done
                    else
                        num_wid = num_wid - multipliers(multid(sortid(j)), 3); 
                        mult_queue(multid(sortid(j))) = mult_queue(multid(sortid(j))) + multipliers(multid(sortid(j)), 3); 
                        
                        if mult_out_ind(multid(sortid(j)), mult_out_ind(multid(sortid(j)),1)+1) == m + out_ind_base*Tm
                            mult_out_num(multid(sortid(j)), mult_out_num(multid(sortid(j)),1)+1) = ...
                                mult_out_num(multid(sortid(j)), mult_out_num(multid(sortid(j)),1)+1) + multipliers(multid(sortid(j)), 3);
                        else
                            mult_out_ind(multid(sortid(j)), 1) = mult_out_ind(multid(sortid(j)), 1) + 1;
                            mult_out_num(multid(sortid(j)), 1) = mult_out_num(multid(sortid(j)), 1) + 1;
                            mult_out_ind(multid(sortid(j)), mult_out_ind(multid(sortid(j)), 1)+1) = m + out_ind_base*Tm;
                            mult_out_num(multid(sortid(j)), mult_out_num(multid(sortid(j)), 1)+1) = multipliers(multid(sortid(j)), 3);
                        end
                        
                    end
                end
            end
            if done_flag
                break;
            else
                dist = dist + 1; 
                if dist > max(multipliers(:,2))
                    break;
                end
                multid = find(multipliers(:,1) == wid-dist);
                multid(multipliers(multid, 3) < dist) = []; 
            end
        end 

    end
           
end