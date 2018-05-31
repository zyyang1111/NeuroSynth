function [mult_use_array, usage] = allocate_mult_use(multipliers, mult_use_array, verbose, blog_assign)

% try to assign the usage of the multiplier
original_mult_num = sum(multipliers(:,3));
if verbose
    fprintf('\nOriginal Unused # Multipliers = %d\n', sum(multipliers(:,3)));
end
for m = 1 : length(mult_use_array)
    if mult_use_array(m).done ~= 0
        continue; 
    end
    mult_use = mult_use_array(m).set;  
    original_num = sum(mult_use(:,2)); % original number of not-assigned weights
    
    save_multipliers = multipliers;
    save_mult_use = mult_use;
    m_done_flag = 1; 
    % after this loop, multipliers saves the unused multipliers
    for i = 1 : size(mult_use, 1)
        wid = mult_use(i,1); 
        num_wid = mult_use(i,2); 
        
        save_wid_num = num_wid; 
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
            % multipliers = save_multipliers; 
            if verbose
                fprintf('Cannot find enough multipliers for %d (within %d assigned %d)\n', wid, save_wid_num, save_wid_num - num_wid);
            end
            m_done_flag = 0;
            if blog_assign
                break;
            end
        end
        
        mult_use(i, 2) = num_wid;
        
    end
    
    
    
    
    if m_done_flag
        mult_use_array(m).set = mult_use; 
        mult_use_array(m).done = 1; 
    else
        if blog_assign
            % if assignment is done for the whole output blog
            mult_use = save_mult_use; 
            multipliers = save_multipliers; 
        end
        mult_use_array(m).set = mult_use; 
        mult_use_array(m).done = 0;
    end
    
    if verbose
        if m_done_flag
            fprintf('Out[%d] Block performs %d multiplications\n', m, original_num - sum(mult_use(:,2)));
            fprintf('Output[%d] is done\t', m);
            fprintf('Rest # Multipliers = %d\n', sum(multipliers(:,3)));
        else 
            fprintf('Out[%d] Block performs %d multiplications\n', m, original_num - sum(mult_use(:,2)));
            fprintf('Output[%d] is NOT done\t', m);
            fprintf('Rest # Multipliers = %d\n', sum(multipliers(:,3)));
        end
    end
                 
end

usage = original_mult_num - sum(multipliers(:,3));
fprintf('Multiplier Usage is %d [ratio = %.4f]\n', usage, (original_mult_num - sum(multipliers(:,3)))/original_mult_num);