function [ cost ] = calc_cost(w_u, idx, str, cluster_size)

if w_u(idx, 1) - w_u(1,1) > cluster_size-1
    cost = inf; 
else
    flag = 1; 
    for i = 2 : idx
        if w_u(i, 3) < w_u(1,3) || w_u(i,4) < w_u(1,4)
            flag = 0;
            break;
        end
    end
    if flag
        if str == 'D'
            cost = w_u(1,3) * idx;
        elseif str == 'A'
            cost = w_u(1,4) * idx;
        else
            cost = inf; 
            fprintf('Unknown flag\n');
        end
    else
        cost = inf; 
    end
end
   
            