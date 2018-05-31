function [mult_queue, mult_out_ind, mult_out_num, add_queue, add_out_ind, add_out_num, tot_energy] = process_mult_new2(multipliers, mult_queue, mult_out_ind, mult_out_num, add_queue, add_out_ind, add_out_num, Tm, tot_energy, out_ind_base)
% this process_mult_new uses the modified out_ind_base

[mult_lib, minW, maxW, len_mult_lib] = generate_mullib(); 
usage = 0; 
for i = 1 : size(multipliers, 1)
    if multipliers(i,2) > 0
        var_flag = 2; 
    else
        var_flag = 1;
    end
    
    num_processed = min(multipliers(i,3), mult_queue(i)); 
    if num_processed == 0
        continue; 
    end
    usage = usage + num_processed; 
    
    % calculate the energy according to the multipliers used 
    energy = sum(max(mult_lib(multipliers(i,1)-minW+1, 2), len_mult_lib(multipliers(i,2)+1, 1)) * num_processed * (mult_lib(multipliers(i,1)-minW+1, 4) + len_mult_lib(multipliers(i,2)+1, 3)));
    energy = energy / 1e6; % [nJ]  
    tot_energy = tot_energy + energy; 
    
    % processed the queue
    mult_queue(i) = mult_queue(i) - num_processed;
    % update mult_out_ind and mult_out_num
    j = 1;
    while j <= mult_out_ind(i,1)
%     for j = 1 : mult_out_ind(i,1)
        if mult_out_num(i, j+1) <= num_processed
            num_processed = num_processed - mult_out_num(i, j+1);
            out_ind = mult_out_ind(i, j+1);
            add_queue(mod(out_ind, 100)) = add_queue(mod(out_ind, 100)) + mult_out_num(i, j+1) * var_flag;    
            if add_out_ind(mod(out_ind, 100), add_out_ind(mod(out_ind, 100),1)+1) == out_ind
                add_out_num(mod(out_ind, 100), add_out_num(mod(out_ind, 100),1)+1) = ...
                    add_out_num(mod(out_ind, 100), add_out_num(mod(out_ind, 100),1)+1) + mult_out_num(i, j+1) * var_flag;
            else
                add_out_ind(mod(out_ind, 100),1) = add_out_ind(mod(out_ind, 100),1) + 1;
                add_out_num(mod(out_ind, 100),1) = add_out_num(mod(out_ind, 100),1) + 1;
                add_out_ind(mod(out_ind, 100), add_out_ind(mod(out_ind, 100),1)+1) = out_ind; 
                add_out_num(mod(out_ind, 100), add_out_num(mod(out_ind, 100),1)+1) = mult_out_num(i, j+1) * var_flag;
            end
            % reset
            for k = j+1 : mult_out_ind(i,1)
                mult_out_ind(i,k) = mult_out_ind(i,k+1);
                mult_out_num(i,k) = mult_out_num(i,k+1);
            end
            mult_out_ind(i, mult_out_ind(i,1)+1) = 0;
            mult_out_num(i, mult_out_num(i,1)+1) = 0;
            mult_out_ind(i,1) = mult_out_ind(i,1) - 1; 
            mult_out_num(i,1) = mult_out_num(i,1) - 1; 
            
        else
            mult_out_num(i,j+1) = mult_out_num(i,j+1) - num_processed; 
            out_ind = mult_out_ind(i,j+1); 
            add_queue(mod(out_ind, 100)) = add_queue(mod(out_ind, 100)) + num_processed * var_flag;
            if add_out_ind(mod(out_ind, 100), add_out_ind(mod(out_ind, 100),1)+1) == out_ind
                add_out_num(mod(out_ind, 100), add_out_num(mod(out_ind, 100),1)+1) = ...
                    add_out_num(mod(out_ind, 100), add_out_num(mod(out_ind, 100),1)+1) + num_processed * var_flag;
            else
                add_out_ind(mod(out_ind, 100),1) = add_out_ind(mod(out_ind, 100),1) + 1;
                add_out_num(mod(out_ind, 100),1) = add_out_num(mod(out_ind, 100),1) + 1;
                add_out_ind(mod(out_ind, 100), add_out_ind(mod(out_ind, 100),1)+1) = out_ind; 
                add_out_num(mod(out_ind, 100), add_out_num(mod(out_ind, 100),1)+1) = num_processed * var_flag;
            end
            
            % reset
            num_processed = 0;
            j = j + 1;
        end
        if num_processed == 0
            break;
        end
    end
end

usage = usage / sum(multipliers(:,3));
% fprintf('Usage is %.3f percent\n', usage*100);