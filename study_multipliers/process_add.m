function [add_queue, add_out_ind, add_out_num, out_done, tot_energy] = process_add(add_queue, add_out_ind, add_out_num, num_add_input, out_done, out_done_num, add_queue_size, add_energy_0, tot_energy)

for i = 1 : size(add_queue, 1)
    num_processed = min(add_queue(i), num_add_input); % num_processed cannot exceed the adder input size
    if num_processed > 0
        num_processed = min(num_processed, add_out_num(i,2)); % num_processed cannot cross multiple outputs
    end
    num_processed = min(num_processed, add_queue_size); % simulate stall due to the limited size of the add-queue
   
    if num_processed == 0
        continue; 
    end
    
    tot_energy = tot_energy + min(fix(num_processed * 1.5), num_add_input) * add_energy_0;
    
    add_queue(i) = add_queue(i) - num_processed; 
    % update
    j = 1;
    while j <= add_out_ind(i,1)
%     for j = 1 : add_out_ind(i,1)
        if add_out_num(i,j+1) <= num_processed
            num_processed = num_processed - add_out_num(i, j+1);
            out_ind = add_out_ind(i, j+1);
            out_done(out_done(:,1) == out_ind, 2) = out_done(out_done(:,1) == out_ind, 2) + add_out_num(i, j+1); 
            % reset
            for k = j+1 : add_out_ind(i,1)
                add_out_ind(i,k) = add_out_ind(i,k+1); 
                add_out_num(i,k) = add_out_num(i,k+1); 
            end
            add_out_ind(i,add_out_ind(i,1)+1) = 0;
            add_out_num(i,add_out_num(i,1)+1) = 0;
            add_out_ind(i,1) = add_out_ind(i,1) - 1;
            add_out_num(i,1) = add_out_num(i,1) - 1;
            
            
        else
            add_out_num(i,j+1) = add_out_num(i,j+1) - num_processed; 
            out_ind = add_out_ind(i,j+1); 
            out_done(out_done(:,1) == out_ind, 2) = out_done(out_done(:,1) == out_ind, 2) + num_processed; 
            
            % reset
            num_processed = 0;
            j = j + 1;
        end
        if num_processed == 0
            break;
        end
    end
end
    