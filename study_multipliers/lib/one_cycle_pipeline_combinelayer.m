function [mult_queue, mult_out_ind, mult_out_num, add_queue, add_out_ind, add_out_num, out_done, cycle_num, tot_energy, stall_cycle, mult_energy, add_energy] = one_cycle_pipeline_combinelayer(multipliers, mult_use, mult_queue, mult_out_ind, mult_out_num, add_queue, add_out_ind, add_out_num, out_ind_base, tm, num_add_input, out_done, out_done_num, cycle_num, mult_queue_size, add_queue_size, add_energy_0, tot_energy, w_map, add_cycle, add_energy, mult_energy)

cycle_num_0 = cycle_num; 
stall_cycle = 0;

[mult_queue, mult_out_ind, mult_out_num] = allocate_mult_use_withmap_combinelayer(multipliers, mult_queue, mult_out_ind, mult_out_num, mult_use, out_ind_base, tm, w_map);
% [mult_queue, mult_out_ind, mult_out_num] = allocate_mult_use_withmap(multipliers, mult_queue, mult_out_ind, mult_out_num, mult_use, out_ind_base, tm, w_map);

while 1
    cycle_num = cycle_num + 1;

    tot_energy_0 = tot_energy;
    [mult_queue, mult_out_ind, mult_out_num, add_queue, add_out_ind, add_out_num, tot_energy] = process_mult_combinelayer(multipliers, mult_queue, mult_out_ind, mult_out_num, add_queue, add_out_ind, add_out_num, tm, tot_energy, out_ind_base);
    mult_energy = mult_energy + tot_energy - tot_energy_0;
    if max(mult_out_ind(:,1))+1 < size(mult_out_ind, 2)
        mult_out_ind(:, max(mult_out_ind(:,1))+2:end) = [];
        mult_out_num(:, max(mult_out_ind(:,1))+2:end) = [];
    end
    
    if max(mult_queue) == 0
        stall_cycle = stall_cycle + 1; 
    end
       
   if mod(cycle_num - cycle_num_0, add_cycle) == 1 || add_cycle == 1
       tot_energy_0 = tot_energy;
       [add_queue, add_out_ind, add_out_num, out_done, tot_energy] = process_add_combinelayer(add_queue, add_out_ind, add_out_num, num_add_input, out_done, out_done_num, add_queue_size, add_energy_0, tot_energy);
       add_energy = add_energy + tot_energy - tot_energy_0;
       if max(add_out_ind(:,1))+1 < size(add_out_ind, 2)
           add_out_ind(:, max(add_out_ind(:,1))+2:end) = [];
           add_out_num(:, max(add_out_ind(:,1))+2:end) = [];
       end
       out_done(out_done(:,2) == out_done_num, :) = [];
   end
   if max(mult_queue - mult_queue_size) < 0 && max(add_queue) < add_queue_size         
       break;
   end

end

% stall_cycle = cycle_num - cycle_num_0 - 1;