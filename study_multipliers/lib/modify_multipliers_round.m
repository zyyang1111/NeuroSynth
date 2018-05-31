function [multipliers] = modify_multipliers_round(multipliers, root_file, suffix)

wr = csvread([root_file 'weight_quant' suffix(7:end) '.csv']);
minW = wr(1);
sizeW = wr(2);
large = wr(3);
quant = wr(4:end);

data = minW : minW+sizeW-1;
quant = [data' , quant];
[mult_lib, minW_lib, maxW, len_mult_lib] = generate_mullib(); 

multipliers_re = multipliers; 
% convert 
for i = 1 : size(multipliers, 1)
    w = multipliers(i,1) : multipliers(i,1) + multipliers(i,2); 
    w_map = quant(w-minW+1, 2); 
    invalid = find(w_map == large); 
    w_map(invalid) = w(invalid); 
    
    base = min(w_map); 
    var = max(w_map) - min(w_map); 
    if base == multipliers(i,1) && var == multipliers(i,2)
        continue; 
    end
    fprintf('Change From (%d, %d) --> (%d, %d)\n', multipliers(i,1), multipliers(i,2), base, var); 
    multipliers(i,1) = base; 
    multipliers(i,2) = var; 
    multipliers(i,3) = multipliers(i,3) + 1;
end
% multipliers(:, 3) = multipliers(:,3) + 1; 

mult_area_re = sum((mult_lib(multipliers_re(:,1)-minW_lib+1, 3) + len_mult_lib(multipliers_re(:,2)+1, 2)) .* multipliers_re(:,3));
mult_area = sum((mult_lib(multipliers(:,1)-minW_lib+1, 3) + len_mult_lib(multipliers(:,2)+1, 2)) .* multipliers(:,3));

fprintf('Original Area: %.4f mm^2\n', mult_area_re / 1e6);
fprintf('New Area: %.4f mm^2\n', mult_area/1e6);
