function [multipliers] = modify_multipliers(multipliers, area_limit)
[mult_lib, minW, maxW, len_mult_lib] = generate_mullib(); 

% modify the multipliers --> combine fixed-weight mult to mult-pair
min_mult = min(multipliers(:,1));
max_mult = max(multipliers(:,1));
mn = [];
for i = min_mult : max_mult
    set = multipliers(multipliers(:,1) == i, :);
    if size(set, 1) == 0
        continue;
    elseif size(set, 1) == 1
        mn = [mn ; set];
    else
        set_n = [set(1,1), max(set(:,2)), sum(set(:,3))];
        mn = [mn ; set_n];
    end
end
multipliers = mn;
% modify 2 
mn = []; 
delete_flag = zeros(size(multipliers, 1), 1);
for i = 1 : size(multipliers, 1)-1
    if delete_flag == 1
        continue;
    end
    base = multipliers(i,1); 
    bias = multipliers(i,2);
    for j = i+1 : size(multipliers, 1)
        if multipliers(j,1)+multipliers(j,2) <= base+bias
            delete_flag(j) = 1; 
        else
            break;
        end
    end
end
for i = 1 : size(multipliers, 1)
    if delete_flag(i) == 0
        mn = [mn ; multipliers(i,:)];
    else
        mn(end, 3) = mn(end, 3) + multipliers(i,3);
    end
end
multipliers = mn;

% modify the number of multipliers to meet the area constraint
area = sum((mult_lib(multipliers(:,1)-minW+1, 3) + len_mult_lib(multipliers(:,2)+1, 2)) .* multipliers(:,3));
A_array = [1000, 800, 500, 300, 100, 50, 20]; 
cnt = 0;
while area > area_limit
    set = find(multipliers(:,3) > A_array(mod(cnt, length(A_array))+1)); 
    multipliers(set, 3) = multipliers(set, 3) - 1; 
    cnt = cnt + 1; 
    area = sum((mult_lib(multipliers(:,1)-minW+1, 3) + len_mult_lib(multipliers(:,2)+1, 2)) .* multipliers(:,3));
    fprintf('Area margin = %d\n', area - area_limit);
end