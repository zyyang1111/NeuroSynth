function [mult_queue, mult_out_ind, mult_out_num, add_queue, add_out_ind, add_out_num, out_done, cycle_num, tot_energy, stall_cycle] = one_cycle_pipeline(multipliers, mult_use, mult_queue, mult_out_ind, mult_out_num, add_queue, add_out_ind, add_out_num, out_ind_base, tm, num_add_input, out_done, out_done_num, cycle_num, mult_queue_size, add_queue_size, add_energy_0, tot_energy)

cycle_num_0 = cycle_num; 

[mult_queue, mult_out_ind, mult_out_num] = allocate_mult_use_new2_new(multipliers, mult_queue, mult_out_ind, mult_out_num, mult_use, out_ind_base, tm);

while 1
    cycle_num = cycle_num + 1;
    
    [mult_queue, mult_out_ind, mult_out_num, add_queue, add_out_ind, add_out_num, tot_energy] = process_mult_new(multipliers, mult_queue, mult_out_ind, mult_out_num, add_queue, add_out_ind, add_out_num, tm, tot_energy);
    if max(mult_out_ind(:,1))+1 < size(mult_out_ind, 2)
        mult_out_ind(:, max(mult_out_ind(:,1))+2:end) = [];
        mult_out_num(:, max(mult_out_ind(:,1))+2:end) = [];
    end
                
                
   [add_queue, add_out_ind, add_out_num, out_done, tot_energy] = process_add(add_queue, add_out_ind, add_out_num, num_add_input, out_done, out_done_num, add_queue_size, add_energy_0, tot_energy);
                
   if max(add_out_ind(:,1))+1 < size(add_out_ind, 2)
       add_out_ind(:, max(add_out_ind(:,1))+2:end) = [];
       add_out_num(:, max(add_out_ind(:,1))+2:end) = [];
   end
   out_done(out_done(:,2) == out_done_num, :) = [];
   if max(mult_queue - mult_queue_size) < 0 && max(add_queue) < add_queue_size
       break;
   end
%    if max(mult_queue) < mult_queue_size && max(add_queue) < add_queue_size
%        break;
%    end
end

stall_cycle = cycle_num - cycle_num_0 - 1;