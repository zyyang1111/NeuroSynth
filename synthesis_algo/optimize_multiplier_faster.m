close all 
clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zhiyuan Yang   May 26 2017
% Try to cluster the multiplications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
root_file = '../../AlexNet_caffe/'; 

Q = 7; 
BIT_SIZE = 16;
area_b = 10071; 
delay_b = 3989;
% load the multiplier lib 
[mult_lib, minW, maxW, len_mult_lib] = generate_mullib(); 


% % Step 1: all the multiplications are separated; find the unique weights
% % and the number of such multiplications 
w2suffix = ''; % this is for rounded w2 file
w_u = generate_combined_wu(root_file, w2suffix);
% w_u = generate_combined_wu_2(root_file, w2suffix, Tm, Tn, Tr);
suffix = ['_comb' w2suffix];
idx = 0;

D = 0; A = 0;
for i = 1 : size(w_u, 1)
    d = mult_lib(w_u(i,1)*2^Q - minW + 1, 2); 
    a = mult_lib(w_u(i,1)*2^Q - minW + 1, 3);
    if d > D
        D = d; 
    end
    A = A + a * w_u(i, 2); 
end
fprintf('======= Initial Case ========\n');
fprintf('Delay = %.4f [ps]\n', D);
fprintf('Area = %.4f [cm^2]\n', A/1e8);

% Step 2: cluster any two multiplications and record the area and delay;
% each cluster is indexed by the type of the combined multiplier; each of
% such multipliers is represented by a two-tuple: (fixed-weight, change-length)
mult = zeros(size(w_u, 1), 5); 
% initialize the mult: store the multiplication information 
% mult(:,1) -- fixed weight 
% mult(:,2) -- length mult (0 if there is no such multiplier)
% mult(:,3) -- num of such multiplier 
% mult(:,4) -- num of cycles of such multiplier to finish all computation
% mult(:,5) -- maximum variance from the baseline (0 if no len-mult)
mult(:,1) = w_u(:,1)*2^Q; mult(:,3) = w_u(:,2); mult(:,4) = 1;
% calculate delay and area using mult:
% delay = max(mult_lib(mult(:,1)-minW+1, 2), len_mult_lib(mult(:,2)+1, 1)) * mult(:,4)
% area = (mult_lib(mult(:,1)-minW+1, 3) + len_mult_lib(mult(:,2)+1, 2)) * mult(:,3)
mult_re = mult;
size_mult = size(mult, 1);
D_array = D; A_array = A; mult_A = [mult ; zeros(1,5)];
full_D_array = D; full_A_array = A; 
for i = 1 : size_mult
    % fprintf('%d\n', i);
    if mult(i,3) < 1
        continue;
    end
    for j = i : size_mult
        if j == i
            if mult(j,3) > 1
                mult_t = mult; 
                mult_t(size_mult+1, 1:2) = mult_t(j,1:2); 
                mult_t(size_mult+1, 3) = fix(mult(j,3) / 2); % faster
                mult_t(j,3) = mult_t(j,3) - 2 * mult_t(size_mult+1, 3); % faster 
                mult_t(size_mult+1, 4) = mult_t(j,4)*2;
                mult_t(size_mult+1, 5) = mult_t(j,5);
                
%                 delay = max(max(mult_lib(mult_t(:,1)-minW+1, 2), len_mult_lib(mult_t(:,2)+1, 1)) .* mult_t(:,4));
                delay = max(max(mult_lib(mult_t(:,1)-minW+1, 2), len_mult_lib(mult_t(:,2)+1, 1))) * max(mult_t(:,4));
                area = sum((mult_lib(mult_t(:,1)-minW+1, 3) + len_mult_lib(mult_t(:,2)+1, 2)) .* mult_t(:,3));
                
                % check the Pareto Optimality 
                discard = 0; 
                remove_array = []; 
                for k = 1 : length(D_array)
                    if D_array(k) <= delay && A_array(k) < area || D_array(k) < delay && A_array(k) <= area
                        discard = 1;
                        break;
                    end
                    if delay <= D_array(k) && area < A_array(k) || delay < D_array(k) && area <= A_array(k)
                        remove_array = cat(1, remove_array, k);
                    end
                end
                
                % add the new results
                if ~discard
                    mult_A = cat(3, mult_A, mult_t); 
                    D_array = cat(1, D_array, delay); 
                    A_array = cat(1, A_array, area); 
                    
                    mult_A(:,:,remove_array) = []; 
                    D_array(remove_array) = []; 
                    A_array(remove_array) = [];
                end

            end
        else
            if mult(j,3) >= 1
                mult_t = mult; 
                m1_l = mult_t(i,1); m1_u = mult_t(i,1) + mult_t(i,5); 
                m2_l = mult_t(j,1); m2_u = mult_t(j,1) + mult_t(j,5); 
                
                ml = min(m1_l, m2_l); 
                mu = max(m1_u, m2_u); 
                if mu == ml
                    mult_t(size_mult+1, 1:2) = [ml, 0]; 
                    mult_t(size_mult+1, 3) = min(mult_t(j,3), mult_t(i,3)); 
                    mult_t(j,3) = mult_t(j,3)-mult_t(size_mult+1, 3); mult_t(i,3) = mult_t(i,3)-mult_t(size_mult+1, 3); 
                    mult_t(size_mult+1, 4) = mult_t(i,4) + mult_t(j,4); 
                    mult_t(size_mult+1, 5) = 0; 
                else
                    mult_t(size_mult+1, 1:2) = [ml, 1+fix(log(mu-ml)/log(2))]; 
                    mult_t(size_mult+1, 3) = min(mult_t(j,3), mult_t(i,3)); 
                    mult_t(j,3) = mult_t(j,3)-mult_t(size_mult+1, 3); mult_t(i,3) = mult_t(i,3)-mult_t(size_mult+1, 3); 
                    mult_t(size_mult+1, 4) = mult_t(i,4) + mult_t(j,4);
                    mult_t(size_mult+1, 5) = mu - ml; 
                end
                
                delay = max(max(mult_lib(mult_t(:,1)-minW+1, 2), len_mult_lib(mult_t(:,2)+1, 1)) .* mult_t(:,4));
                area = sum((mult_lib(mult_t(:,1)-minW+1, 3) + len_mult_lib(mult_t(:,2)+1, 2)) .* mult_t(:,3));
                
                full_D_array = cat(1, full_D_array, delay);
                full_A_array = cat(1, full_A_array, area);
                
                % check the Pareto Optimality 
                discard = 0; 
                remove_array = []; 
                for k = 1 : length(D_array)
                    if D_array(k) <= delay && A_array(k) < area || D_array(k) < delay && A_array(k) <= area
                        discard = 1;
                        break;
                    end
                    if delay <= D_array(k) && area < A_array(k) || delay < D_array(k) && area <= A_array(k)
                        remove_array = cat(1, remove_array, k);
                    end
                end
                
                % add the new results
                if ~discard
                    mult_A = cat(3, mult_A, mult_t); 
                    D_array = cat(1, D_array, delay); 
                    A_array = cat(1, A_array, area); 
                    
                    mult_A(:,:,remove_array) = []; 
                    D_array(remove_array) = []; 
                    A_array(remove_array) = [];
                end
            end
        end
    end
end

% % check som consistency 
% fprintf('========== Step 2 ===========\n');
% fprintf('Original Number of Multipliers is %d\n', sum(mult_re(:,3))); 
% fprintf('Original Number Multiplications is %d\n', sum(mult_re(:,3) .* mult_re(:,4)));
% for l = 1 : size(mult_A, 3)
%     fprintf('%d: NMer: %d; NM: %d; Delay = %.2f; Area = %.2f\n', l, sum(mult_A(:,3, l)), sum(mult_A(:,3, l) .* mult_A(:,4, l)), D_array(l), A_array(l));
% end

fprintf('========== Step 2 ===========\n');
[minA, minA_ind] = min(A_array); 
fprintf('min Area = %d; Delay = %d\n', minA, D_array(minA_ind));

% figure
% plot(A_array, D_array, 'b*');
% hold on 

% Step 3: Now we have a set of pareto multiplication: mult_A: [a x b x c];
% we will check the 2-D array one by one and stores the Pareto Optimal
% results for this level 
% The methods is similar to the Step 2 except for another dimension is
% considered
step = 2; 


while 1
    step = step + 1; 
mult_A = cat(1, mult_A, zeros(1, 5, size(mult_A, 3)));
mult_A_re = mult_A; 
size_mult = size(mult_A, 1) - 1; 
for l = 1 : size(mult_A_re, 3)
    mult = mult_A_re(1:end-1, :, l); 
    for i = 1 : size_mult
        % fprintf('%d\n', i);
        if mult(i,3) < 1
            continue;
        end
        for j = i : size_mult
            if j == i
                if mult(j,3) > 1
                    mult_t = mult; 
                    mult_t(size_mult+1, 1:2) = mult_t(j,1:2); 
                    mult_t(size_mult+1, 3) = fix(mult(j,3) / 2); % faster
                    mult_t(j,3) = mult_t(j,3) - 2 * mult_t(size_mult+1, 3); % faster 
                    mult_t(size_mult+1, 4) = mult_t(j,4)*2;
                    mult_t(size_mult+1, 5) = mult_t(j,5);

                    % check if we can combine two multipliers 
                    for k = size_mult : -1 : 1
                        if mult_t(k,1) == mult_t(size_mult+1,1) && mult_t(k,2) == mult_t(size_mult+1,2) && mult_t(k,4) == mult_t(size_mult+1,4) && mult_t(k,5) == mult_t(size_mult+1,5)
                            mult_t(k,3) = mult_t(k,3) + mult_t(size_mult+1,3); 
                            mult_t(size_mult+1, :) = 0;
                            break;
                        end
                    end
                    % check the non-used multipliers
                    remove_array = [];
                    for k = 1 : size_mult+1
                        if mult_t(k,3) < 1
                            remove_array = cat(1, remove_array, k); 
                        end
                    end
                    if ~isempty(remove_array)
                        mult_t(remove_array, :) = [];
                        mult_t = cat(1, mult_t, zeros(length(remove_array), 5));
                    end
                    
                    if sum(mult_t(:,3).*mult_t(:,4)) ~= 105415200
                        call_stop = 1;
                    end
                
                    
                    delay = max(max(mult_lib(mult_t(:,1)-minW+1, 2), len_mult_lib(mult_t(:,2)+1, 1))) * max(mult_t(:,4));
                    area = sum((mult_lib(mult_t(:,1)-minW+1, 3) + len_mult_lib(mult_t(:,2)+1, 2)) .* mult_t(:,3));
                
                    % check the Pareto Optimality 
                    discard = 0; 
                    remove_array = []; 
                    for k = 1 : length(D_array)
                        if D_array(k) <= delay && A_array(k) <= area || D_array(k) <= delay && A_array(k) <= area
                            discard = 1;
                            break;
                        end
                        if delay <= D_array(k) && area < A_array(k) || delay < D_array(k) && area <= A_array(k)
                            remove_array = cat(1, remove_array, k);
                        end
                    end
                
                    % add the new results
                    if ~discard
                        mult_A = cat(3, mult_A, mult_t); 
                        D_array = cat(1, D_array, delay); 
                        A_array = cat(1, A_array, area); 
                    
                        mult_A(:,:,remove_array) = []; 
                        D_array(remove_array) = []; 
                        A_array(remove_array) = [];
                    end
                end
            else
                if mult(j,3) >= 1
                    mult_t = mult; 
                    m1_l = mult_t(i,1); m1_u = mult_t(i,1) + mult_t(i,5); 
                    m2_l = mult_t(j,1); m2_u = mult_t(j,1) + mult_t(j,5); 
                
                    ml = min(m1_l, m2_l); 
                    mu = max(m1_u, m2_u); 
                    if mu-ml <= 5
                        if mu == ml
                            mult_t(size_mult+1, 1:2) = [ml, 0]; 
                            mult_t(size_mult+1, 3) = min(mult_t(j,3), mult_t(i,3)); 
                            mult_t(j,3) = mult_t(j,3)-mult_t(size_mult+1, 3); mult_t(i,3) = mult_t(i,3)-mult_t(size_mult+1, 3); 
                            mult_t(size_mult+1, 4) = mult_t(i,4) + mult_t(j,4); 
                            mult_t(size_mult+1, 5) = 0; 
                        else
                            mult_t(size_mult+1, 1:2) = [ml, 1+fix(log(mu-ml)/log(2))]; 
                            mult_t(size_mult+1, 3) = min(mult_t(j,3), mult_t(i,3)); 
                            mult_t(j,3) = mult_t(j,3)-mult_t(size_mult+1, 3); mult_t(i,3) = mult_t(i,3)-mult_t(size_mult+1, 3); 
                            mult_t(size_mult+1, 4) = mult_t(i,4) + mult_t(j,4);
                            mult_t(size_mult+1, 5) = mu - ml; 
                        end
                    
                        % check if we can combine two multipliers 
                        for k = size_mult : -1 : 1
                            if mult_t(k,1) == mult_t(size_mult+1,1) && mult_t(k,2) == mult_t(size_mult+1,2) && mult_t(k,4) == mult_t(size_mult+1,4) && mult_t(k,5) == mult_t(size_mult+1,5)
                                mult_t(k,3) = mult_t(k,3) + mult_t(size_mult+1,3); 
                                mult_t(size_mult+1, :) = 0;
                                break;
                            end
                        end
                        % check the non-used multipliers
                        remove_array = [];
                        for k = 1 : size_mult+1
                            if mult_t(k,3) < 1
                                remove_array = cat(1, remove_array, k); 
                            end
                        end
                        if ~isempty(remove_array)
                            mult_t(remove_array, :) = [];
                            mult_t = cat(1, mult_t, zeros(length(remove_array), 5));
                        end
                      
                        if sum(mult_t(:,3).*mult_t(:,4)) ~= 105415200
                            call_stop = 1;
                        end
                    
                        delay = max(max(mult_lib(mult_t(:,1)-minW+1, 2), len_mult_lib(mult_t(:,2)+1, 1))) * max(mult_t(:,4));
                        area = sum((mult_lib(mult_t(:,1)-minW+1, 3) + len_mult_lib(mult_t(:,2)+1, 2)) .* mult_t(:,3));
                
                        % check the Pareto Optimality 
                        discard = 0; 
                        remove_array = []; 
                        for k = 1 : length(D_array)
                            if D_array(k) <= delay && A_array(k) <= area || D_array(k) <= delay && A_array(k) <= area
                                discard = 1;
                                break;
                            end
                            if delay <= D_array(k) && area < A_array(k) || delay < D_array(k) && area <= A_array(k)
                                remove_array = cat(1, remove_array, k);
                            end
                        end
                
                        % add the new results
                        if ~discard
                            mult_A = cat(3, mult_A, mult_t); 
                            D_array = cat(1, D_array, delay); 
                            A_array = cat(1, A_array, area); 
                    
                            mult_A(:,:,remove_array) = []; 
                            D_array(remove_array) = []; 
                            A_array(remove_array) = [];
                        end
                    end
                end
            end
        end
    end
end

% delete the pareto optimal results with large area
% truncate the pareto optimal size to speed up the simulation 
minA = min(A_array);
fprintf('minA = %.6f; maxA = %.2f\n', minA/1e8, max(A_array)/1e8);
large_A = find(A_array > minA * 3); % find the case when the area is too large
A_array(large_A) = [];
D_array(large_A) = [];
mult_A(:,:,large_A) = [];


% delete the zero lines in mult_A
mult_A_sum = sum(sum(abs(mult_A), 3), 2); 
zero_mult_A_sum = find(mult_A_sum == 0); 
if ~isempty(zero_mult_A_sum)
    mult_A(zero_mult_A_sum, :, :) = [];
end

% % check consistency 
fprintf('========== Step %d ===========\n', step);
fprintf('Remove %d pareto results\n', length(large_A));
fprintf('mult_A size is %d; Pareto Size is %d\n', size(mult_A, 1), size(mult_A, 3));
[minA, minA_ind] = min(A_array); 
fprintf('min Area = %.6f [cm^2]; Delay = %.2f [ps]\n', minA/1e8, D_array(minA_ind));


% Set the saving points for different area constraints
if minA/1e8 <= 0.0001
    save([root_file num2str(idx) 'result_00001' suffix '.mat']);
    break;
elseif minA/1e8 <= 0.0005
    save([root_file num2str(idx) 'result_00005' suffix '.mat']);
elseif minA/1e8 <= 0.001
    save([root_file num2str(idx) 'result_0001' suffix '.mat']);
elseif minA/1e8 <= 0.003
    save([root_file num2str(idx) 'result_0003' suffix '.mat']);
elseif minA/1e8 <= 0.006
    save([root_file num2str(idx) 'result_0006' suffix '.mat']);
elseif minA/1e8 <= 0.008
    save([root_file num2str(idx) 'result_0008' suffix '.mat']);
elseif minA/1e8 <= 0.01
    save([root_file num2str(idx) 'result_001' suffix '.mat']);
elseif minA/1e8 <= 0.012
    save([root_file num2str(idx) 'result_0012' suffix '.mat']);
elseif minA/1e8 <= 0.015
    save([root_file num2str(idx) 'result_0015' suffix '.mat']);
elseif minA/1e8 <= 0.018
    save([root_file num2str(idx) 'result_0018' suffix '.mat']);
elseif minA/1e8 <= 0.021
    save([root_file num2str(idx) 'result_0021' suffix '.mat']);
elseif minA/1e8 <= 0.025
    save([root_file num2str(idx) 'result_0025' suffix '.mat']);
elseif minA/1e8 <= 0.03
    save([root_file num2str(idx) 'result_003' suffix '.mat']);
elseif minA/1e8 <= 0.035
    save([root_file num2str(idx) 'result_0035' suffix '.mat']);
end

end
% plot(A_array, D_array, 'r*');    
