close all 
clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zhiyuan Yang  June 8 2017 
% Check whether the current multiplier lib is enough
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
root_file = '../../Lenet/'; 
addpath './num_convert/';

Q = 7; 
BIT_SIZE = 16;
area_b = 10071; 
delay_b = 3989;
% load the multiplier lib 
[mult_lib, minW, maxW, len_mult_lib] = generate_mullib(); 

[fig_size, file_array] = read_NetInfo(root_file);

weights = [];
num_weights = zeros(length(fig_size), 1);
num_syn = zeros(length(fig_size), 1);
for idx = 1 : length(fig_size)
    file = [root_file file_array{idx} '_w2.csv'];
    w = csvread(file);
    w_u = unique(w); 
    num_weights(idx) = length(w_u);
    num_syn(idx) = length(w(:))*fig_size(idx)^2;
    
    weights = union(weights, w_u);
    fprintf('Layer %d: num_weights = %d; num_syn = %d\n', idx, num_weights(idx), num_syn(idx));
end

fprintf('\nTotal number of different weights = %d\n', length(weights));

fprintf('============= check lib ===========\n');
flag_exceed = 0; 
flag_miss = 0; 
for i = 1 : length(weights)
    w = weights(i)*2^Q;
    if w < minW || w > maxW 
        flag_exceed = 1;
        fprintf('%d\n', w);
    else
        if mult_lib(w-minW+1, 3) == 0
            flag_miss = 1;
            fprintf('%d\n', w);
        end
    end
end
if ~flag_exceed && ~flag_miss 
    fprintf('All the weights are in the library!\n');
end
