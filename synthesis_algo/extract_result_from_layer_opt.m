close all 
clear
clc

% Zhiyuan Yang   May 27 2017
% get the specific result from the results optimized for each layer
root_file = '../../GenderNet/'; 
% load the multiplier lib 
[mult_lib, minW, maxW, len_mult_lib] = generate_mullib(); 

Alimit = [0.003 0.006 0.01 0.012];
Pareto = [0.003 0.006 0.01 0.012];


w2suffix = '_round_2';
suffix = ['_comb' w2suffix];

for iter = 1 : length(Alimit)
    area_limit = Alimit(iter); 
    Pareto_set = num2str(Pareto(iter)); 
    Pareto_set(Pareto_set == '.') = []; 
    save_set = num2str(area_limit); 
    save_set(save_set == '.') = [];


% Pareto_set = '001'; 
% area_limit = 0.01; 
% save_set = '001';
save_file = [root_file 'separate_results/P' Pareto_set '_A' save_set suffix '/'];
mkdir(save_file);

% load the network information 
% [fig_size, file_array] = determine_NetInfo('AlexNet_caffe');
[fig_size, file_array] = read_NetInfo(root_file);

save_delay = zeros(length(fig_size) , 1);
save_clk = zeros(length(fig_size) , 1);
save_cyc = zeros(length(fig_size) , 1);

for i = 0 : 0
% for i = 1 : length(fig_size)
    file = [root_file num2str(i) 'result_' Pareto_set suffix '.mat'];
    load(file, 'D_array', 'A_array', 'mult_A'); 
    [sort_A, sort_ID] = sort(A_array); 
    sort_D = D_array(sort_ID); 
    
    flag = 0;
    for j = 1 : length(sort_A)
        if sort_A(j) > area_limit * 1e8 && j > 1
            idx = j - 1;
            flag = 1; 
            break;
        end
    end
    if ~flag 
        idx = length(sort_A);
    end
    
    area = sort_A(idx) / 1e8;
    delay = sort_D(idx) / 1e6;
    mult_input = mult_A(:,:,sort_ID(idx));
    
    clk = max(max(mult_lib(mult_input(:,1)-minW+1, 2), len_mult_lib(mult_input(:,2)+1, 1)));
    cyc = max(mult_input(:,4));
    delay_n = clk * cyc;
    
    if i == 0
        jdx = 1;
    else
        jdx = i;
    end
    save([save_file 'mult_layer' num2str(i) '_' save_set '.mat'], 'mult_input');
    save_delay(jdx) = delay;
    save_clk(jdx) = clk;
    save_cyc(jdx) = cyc;
    fprintf('Layer %d: Area = %.6f [cm^2];  Delay = %.6f [us]\n', i, area, delay);
    fprintf('Total Multiplications is %d\n', sum(mult_input(:,3).*mult_input(:,4)));
    fprintf('Clk = %.2f [ns], #cycles = %d, delay = %.6f [us]\n\n', clk/1e3, cyc, delay_n/1e6);
end

fprintf('\n');
end