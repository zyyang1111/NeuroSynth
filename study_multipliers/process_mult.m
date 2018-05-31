function [mult_queue, add_queue, add_out_ind, add_out_num] = process_mult(multipliers, mult_queue, add_queue, add_out_ind, add_out_num, Tm)

for i = 1 : size(multipliers, 1)
    num_processed = min(multipliers(i,3), mult_queue(i).queue); 
    if num_processed == 0
        continue; 
    end
    
    % processed the queue
    mult_queue(i).queue = mult_queue(i).queue - num_processed;
    % update mult_out_ind and mult_out_num
    j = 1;
    while j <= mult_queue(i).out_ind(1)
%     for j = 1 : mult_out_ind(i,1)
        if mult_queue(i).out_num(j+1) <= num_processed
            num_processed = num_processed - mult_queue(i).out_num(j+1);
            out_ind = mult_queue(i).out_ind(j+1);
            add_queue(mod(out_ind-1, Tm)+1) = add_queue(mod(out_ind-1, Tm)+1) + mult_queue(i).out_num(j+1);
            if add_out_ind(mod(out_ind-1, Tm)+1, add_out_ind(mod(out_ind-1, Tm)+1,1)+1) == out_ind
                add_out_num(mod(out_ind-1, Tm)+1, add_out_num(mod(out_ind-1, Tm)+1,1)+1) = ...
                    add_out_num(mod(out_ind-1, Tm)+1, add_out_num(mod(out_ind-1, Tm)+1,1)+1) + mult_queue(i).out_num(j+1);
            else
                add_out_ind(mod(out_ind-1, Tm)+1,1) = add_out_ind(mod(out_ind-1, Tm)+1,1) + 1;
                add_out_num(mod(out_ind-1, Tm)+1,1) = add_out_num(mod(out_ind-1, Tm)+1,1) + 1;
                add_out_ind(mod(out_ind-1, Tm)+1, add_out_ind(mod(out_ind-1, Tm)+1,1)+1) = out_ind; 
                add_out_num(mod(out_ind-1, Tm)+1, add_out_num(mod(out_ind-1, Tm)+1,1)+1) = mult_queue(i).out_num(j+1);
            end
            % reset
            for k = j+1 : mult_queue(i).out_ind(1)
                mult_queue(i).out_ind(k) = mult_queue(i).out_ind(k+1);
                mult_queue(i).out_num(k) = mult_queue(i).out_num(k+1);
            end
            mult_queue(i).out_ind(mult_queue(i).out_ind(1)+1) = 0;
            mult_queue(i).out_num(mult_queue(i).out_num(1)+1) = 0;
            mult_queue(i).out_ind(1) = mult_queue(i).out_ind(1) - 1; 
            mult_queue(i).out_num(1) = mult_queue(i).out_num(1) - 1; 
            
        else
            mult_queue(i).out_num(j+1) = mult_queue(i).out_num(j+1) - num_processed; 
            out_ind = mult_queue(i).out_ind(j+1); 
            add_queue(mod(out_ind-1, Tm)+1) = add_queue(mod(out_ind-1, Tm)+1) + num_processed;
            if add_out_ind(mod(out_ind-1, Tm)+1, add_out_ind(mod(out_ind-1, Tm)+1,1)+1) == out_ind
                add_out_num(mod(out_ind-1, Tm)+1, add_out_num(mod(out_ind-1, Tm)+1,1)+1) = ...
                    add_out_num(mod(out_ind-1, Tm)+1, add_out_num(mod(out_ind-1, Tm)+1,1)+1) + num_processed;
            else
                add_out_ind(mod(out_ind-1, Tm)+1,1) = add_out_ind(mod(out_ind-1, Tm)+1,1) + 1;
                add_out_num(mod(out_ind-1, Tm)+1,1) = add_out_num(mod(out_ind-1, Tm)+1,1) + 1;
                add_out_ind(mod(out_ind-1, Tm)+1, add_out_ind(mod(out_ind-1, Tm)+1,1)+1) = out_ind; 
                add_out_num(mod(out_ind-1, Tm)+1, add_out_num(mod(out_ind-1, Tm)+1,1)+1) = num_processed;
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