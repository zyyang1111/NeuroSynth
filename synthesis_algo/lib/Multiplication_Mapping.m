function [clk_new, time, totE, new_power] = Multiplication_Mapping(root_root_file, file_str, str, suffix, w2suffix, sample_layer)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zhiyuan Yang June 1 2017
% Assign the layer1 results to other layers and calculate the delay 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% root_root_file = '../../VGG_19/'; 
Q = 7; 
BIT_SIZE = 16;
area_b = 10071; 
delay_b = 3989;
delay_add = 1702;
% load the multiplier lib 
[mult_lib, minW, maxW, len_mult_lib] = generate_mullib(); 

% [fig_size, file_array] = determine_NetInfo('AlexNet_caffe');
[fig_size, file_array] = read_NetInfo(root_root_file);

num_layer = length(fig_size);

% str = '001'; 
% suffix = '_comb';
% sample_layer = 0;

root_file = [root_root_file 'separate_results/P' file_str '_A' str suffix '/']; 
 
load([root_file 'mult_layer' num2str(sample_layer) '_' str '.mat']);
mkdir([root_file 'total_results/']);

% Step 1: sort the multiplier according to the weight 
[~, ID] = sort(mult_input(:,1)); 
mult_input = mult_input(ID, :); 

% Step 2: calculate the resources of the multipliers 
% if mult_layer1(:,1) and mult_layer1(:,2) are the same, they are the same multiplier 
% sum up mult_layer1(:,3)
multipliers = mult_input(1, 1:3);
for i = 2 : size(mult_input, 1)
    flag = 0;
    for j = 1 : size(multipliers, 1)
        if mult_input(i,1) == multipliers(j,1) && mult_input(i,2) == multipliers(j, 2)
            multipliers(j,3) = multipliers(j,3) + mult_input(i,3);
            flag = 1; 
            break;
        end
    end
    if ~flag
        multipliers = cat(1, multipliers, mult_input(i, 1:3));
    end
end
delay_re = max(max(mult_lib(mult_input(:,1)-minW+1, 2), len_mult_lib(mult_input(:,2)+1, 1))) * max(mult_input(:,4));
area_re = sum((mult_lib(multipliers(:,1)-minW+1, 3) + len_mult_lib(multipliers(:,2)+1, 2)) .* multipliers(:,3));
fprintf('Original Area is %.5f [cm^2]\n', area_re/1e8);
fprintf('Layer 1 Delay is %.5f [us]\n\n', delay_re/1e6);

area_save = zeros(num_layer, 1);
delay_save = zeros(num_layer, 1);
clk_save = zeros(num_layer, 1);
cyc_save = zeros(num_layer, 1);
ene_save = zeros(num_layer, 1);



for idx = 1 : length(fig_size)
% for idx = 1 : num_layer
    num_save = 0;
    mult_sum_save = 0;
    
    file = [root_root_file file_array{idx} '_w2' w2suffix '.csv'];
    w_2D = csvread(file);
    w = w_2D(:);
    
    w_u = unique(w); 
    for i = 1 : size(w_u, 1)
        w_u(i,2) = length(find(w == w_u(i,1))) * fig_size(idx)^2;
    end
    % convert the weight to integers 
    w_u(:,1) = w_u(:,1) * 2^Q;
    
    % build the connection matrix 
    con_mat = zeros(size(multipliers, 1), size(w_u, 1)); 
    for i = 1 : size(multipliers, 1)
        for j = 1 : size(w_u, 1)
            if w_u(j) >= multipliers(i,1) && w_u(j) <= multipliers(i,1)+multipliers(i,2)
                con_mat(i,j) = 1;
            end
        end
    end
    
    % find the set of multiplications that cannot map to the current
    % multipliers 
    new_mult = sum(con_mat, 1); 
    new_mult = find(new_mult == 0);
    
    mult_match = zeros(size(multipliers, 1), size(w_u, 1));
    
    A = con_mat * con_mat';
    A(A > 0) = 1; 
    mult_set = zeros(size(multipliers, 1), 1);
    mult_assign = zeros(size(w_u, 1), 1);
    mult = [];
    
    wu_set = w_u(:,1)*0;
    for i = 1 : size(multipliers, 1)
        if ~mult_set(i)
            len_cl = 0;
            cl = find(A(:,i) > 0); 
            while length(cl) ~= len_cl
                len_cl = length(cl);
                ncl = [];
                for j = 1 : length(cl)
                    ncl = cat(1, ncl, find(A(:,cl(j)) > 0)); 
                end
                cl = union(cl, ncl);
            end
            
            for j = 1 : length(cl)
                mult_set(cl(j)) = 1; 
            end
            
            subset = con_mat(cl, :);
            n_mult = zeros(length(cl), 1); % storing the number of multiplications 
            for k = 1 : size(w_u, 1)
                if sum(subset(:,k)) == 0
                    continue;
                elseif sum(subset(:,k)) == 1
                    n_mult(subset(:,k) == 1) = n_mult(subset(:,k) == 1) + w_u(k,2); 
                    mult_match(cl(subset(:,k) == 1), k) = w_u(k,2);
                    mult_sum_save = mult_sum_save + w_u(k,2);
                end
            end
            for k = 1 : size(w_u, 1)
                if sum(subset(:,k)) > 1
                    mult_sum_save = mult_sum_save + w_u(k,2);
                    mult_sel = find(subset(:,k) == 1); 
                    a = multipliers(cl(mult_sel), 3); 
                    b = n_mult(mult_sel); 
                    f = zeros(length(mult_sel)+1, 1); 
                    f(end) = 1; 
                    aa = zeros(length(mult_sel), length(mult_sel)+1); 
                    for j = 1 : length(mult_sel)
                        aa(j,j) = 1; 
                        aa(j,length(mult_sel)+1) = -a(j);
                    end
                    aeq = ones(1,length(mult_sel)+1); aeq(end)=0;
                    beq = w_u(k, 2); 
                    intcon = 1:length(mult_sel); 
                    options = optimoptions('intlinprog','Display','off');
                    x = intlinprog(f, intcon, aa, -b, aeq, beq, [], [], options); 
                    
                    mult_match(cl(mult_sel), k) = round(x(1:end-1));
                    n_mult(mult_sel) = n_mult(mult_sel) + x(1:end-1);
                end
            end
            
            n_mult = round(n_mult);
            
            yzy = sum(subset, 1); 
            nop = sum(w_u(yzy > 0, 2));
            wu_set(yzy>0) = wu_set(yzy>0) + 1;
            if nop ~= sum(n_mult)
                fprintf('error\n');
            end
            
            
            
            num_save = num_save + sum(n_mult);
%             fprintf('sum_1 = %d; sum_2 = %d\n', mult_sum_save, num_save);
            
            for j = 1 : length(cl)
                mper = cl(j); 
                [NMer, NM] = divide_mult(multipliers(mper, 3), n_mult(j));
                for k = 1 : length(NM)
                    mult_t = [multipliers(mper, 1), multipliers(mper, 2), NMer(k), NM(k), 1]; 
                    mult = cat(1, mult, mult_t); 
                end
            end
        end
    end
    
    % calculate the delay of the current usage 
    d0 = max(max(mult_lib(multipliers(:,1)-minW+1, 2), len_mult_lib(multipliers(:,2)+1, 1))) * max(mult(:,4));
    maxd0 = max(max(mult_lib(multipliers(:,1)-minW+1, 2), len_mult_lib(multipliers(:,2)+1, 1))); 
    maxmult0 = max(mult(:,4));
    % add new multipliers 
    for ii = 1 : length(new_mult)
        i = new_mult(ii);
        n0 = w_u(i,2); 
        [NMer, NM] = assign_newmult(maxmult0, n0); 
        for j = 1 : length(NM)
            mult_t = [w_u(i,1), 0, NMer(j), NM(j), 0];
            mult = cat(1, mult, mult_t);
        end
        multipliers_t = [w_u(i,1), 0, sum(NMer)]; 
        multipliers = cat(1, multipliers, multipliers_t);
        num_save = num_save + n0;
        mult_match(size(mult_match,1)+1, i) = w_u(i,2);
    end
    
    delay = max(max(mult_lib(multipliers(:,1)-minW+1, 2), len_mult_lib(multipliers(:,2)+1, 1))) * max(mult(:,4));
    
    energy = sum(max(mult_lib(mult(:,1)-minW+1, 2), len_mult_lib(mult(:,2)+1, 1)).*mult(:,4).*mult(:,3).*(mult_lib(mult(:,1)-minW+1, 4) + len_mult_lib(mult(:,2)+1, 3)));
    energy = energy / 1e6; % [nJ]  
    
    area = sum((mult_lib(multipliers(:,1)-minW+1, 3) + len_mult_lib(multipliers(:,2)+1, 2)) .* multipliers(:,3));
    clk = max(max(mult_lib(multipliers(:,1)-minW+1, 2), len_mult_lib(multipliers(:,2)+1, 1)));
    cyc = max(mult(:,4));
    
    
    area_save(idx) = area;
    delay_save(idx) = delay; 
    clk_save(idx) = clk;
    cyc_save(idx) = cyc;
    ene_save(idx) = energy;

    fprintf('Layer %d num_save = %d, Total multiplications is %d\n', idx, sum(w_u(:,2)), sum(mult(:,3).*mult(:,4)));
    fprintf('Layer %d Area is %.5f [cm^2]\n', idx, area/1e8);
    fprintf('Layer %d Delay is %.5f [us]\n', idx, delay/1e6);
    fprintf('Layer %d Energy is %d [nJ]\n', idx, energy);
    fprintf('Layer %d Clk is %d [ps], cyc is %d\n\n', idx, clk, cyc);
    
    save([root_file 'total_results/layer_' num2str(idx) '_mult.mat'], 'mult', 'w_2D', 'w_u', 'mult_match');
    
end



delay_save = delay_save / 1e6;
delay_new = max(delay_add, max(clk_save)) * cyc_save / 1e6;
clk_new = max(delay_add, max(clk_save));


power = sum((mult_lib(multipliers(:,1)-minW+1, 4).*mult_lib(multipliers(:,1)-minW+1, 2)/2040 + ...
    len_mult_lib(multipliers(:,2)+1, 3).*len_mult_lib(multipliers(:,2)+1, 3)/2040) .* multipliers(:,3));


sum(multipliers(:,3));


time = sum(delay_new);
totE = sum(ene_save);

new_power = totE * 1e-9 / (time * 1e-6) * 1e3;

save([root_file 'total_results/multipliers.mat'], 'multipliers', 'clk_new');
