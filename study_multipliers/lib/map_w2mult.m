function [w_u] = map_w2mult(w, multipliers)

w_u = unique(w); 

for i = 1 : size(w_u, 1)
    wid = w_u(i);
    
    if wid == 54
        fprintf('check\n');
    end
    
    dist = 0; 
    multid = find(multipliers(:,1) == wid - dist);
    multid(multipliers(multid, 3) < dist) = []; 
    done_flag = 0;
    while ~done_flag
        if ~isempty(multid)
            [~, sortid] = sort(multipliers(multid, 2)); 
            w_u(i,2) = multid(sortid(1));
            done_flag = 1; 
        else
            dist = dist + 1;
            if dist > max(multipliers(:,2))
                fprintf('Cannot find the multiplier for %d\n', wid); 
                break;
            end
            multid = find(multipliers(:,1) == wid-dist);
            multid(multipliers(multid, 2) < dist) = [];
        end
    end
end