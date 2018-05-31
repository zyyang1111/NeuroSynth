function [ cost , sel , array ] = calc_cost_new(w_u, idx, str, cluster_size)
% This cost calculation function will find the centroid within the
% clustered weight rather than using the root value as in "calc_cost_new"

if w_u(idx, 1) - w_u(1,1) > cluster_size-1
    cost = inf; 
    sel = [];
    array = [];
else
    flag = 0; 
    D = w_u(1,3); A = w_u(1,4);sel = 1;
    for i = 1 : idx
        if w_u(i, 3) <= D && w_u(i,4) <= A
            flag = 1;
            D = w_u(i,3); A = w_u(i,4); sel = i;
        end
    end
    if flag
        if str == 'D'
            cost = D * idx;
            array = 1 : idx;
        elseif str == 'A'
            cost = A * idx;
            array = 1 : idx;
        else
            cost = inf; 
            sel = [];
            array = [];
            fprintf('Unknown flag\n');
        end
    else
        cost = inf; 
    end
end