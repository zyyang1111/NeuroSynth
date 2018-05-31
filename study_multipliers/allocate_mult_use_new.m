function [mult_use_array, total_done, mult_left] = allocate_mult_use_new(multipliers, mult_use_array)

total_done = 1; 
mult_left = 0;
% try to assign the usage of the multiplier
original_mult_num = sum(multipliers(:,3));

for m = 1 : length(mult_use_array)
    if mult_use_array(m).done ~= 0
        continue; 
    end
    mult_use = mult_use_array(m).set;  
    
    m_done_flag = 1; 
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
                        multipliers(multid(sortid(j)), 3) = multipliers(multid(sortid(j)), 3) - num_wid; 
                        num_wid = 0;
                        done_flag = 1; 
                        break;
                        % done
                    else
                        num_wid = num_wid - multipliers(multid(sortid(j)), 3); 
                        multipliers(multid(sortid(j)), 3) = 0;
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
        
        if ~done_flag 
            m_done_flag = 0;
        end
        
        mult_use(i, 2) = num_wid;
        
    end
    
    if m_done_flag
        mult_use_array(m).set = mult_use; 
        mult_left = mult_left + sum(mult_use(:,2));
        mult_use_array(m).done = 1; 
    else
        mult_use_array(m).set = mult_use; 
        mult_left = mult_left + sum(mult_use(:,2));
        mult_use_array(m).done = 0;
        total_done = 0;
    end
                 
end

usage = original_mult_num - sum(multipliers(:,3));
fprintf('Multiplier Usage is %d [ratio = %.4f]\n', usage, (original_mult_num - sum(multipliers(:,3)))/original_mult_num);