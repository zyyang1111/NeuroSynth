function [add_queue, add_out_ind, add_out_num, out_done, tot_energy] = empty_outbuffer(add_queue, add_out_ind, add_out_num, num_add_input, out_done, out_done_num, add_queue_size, add_energy_0, tot_energy)

while ~isempty(out_done) && max(add_queue) > 0
    [add_queue, add_out_ind, add_out_num, out_done, tot_energy] = process_add(add_queue, add_out_ind, add_out_num, num_add_input, out_done, out_done_num, add_queue_size, add_energy_0, tot_energy);
                
   if max(add_out_ind(:,1))+1 < size(add_out_ind, 2)
       add_out_ind(:, max(add_out_ind(:,1))+2:end) = [];
       add_out_num(:, max(add_out_ind(:,1))+2:end) = [];
   end
   out_done(out_done(:,2) == out_done_num, :) = [];
end