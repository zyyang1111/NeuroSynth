close all 
clear
clc

addpath './parse_weight/';

root_file = '../AlexNet_caffe/'; 
fr = fopen([root_file 'alexnet_weight.txt'], 'r');

Q = 7;

layer_cnt = 0; 
while ~feof(fr)
    line = fgetl(fr); 
    if line(1) == 'R'
        layer_cnt = layer_cnt + 1; 
        num_line = fgetl(fr); 
        digit_idx = isstrprop(num_line, 'digit');
        num = str2num(num_line(digit_idx == 1));
        data = zeros(num, 1); 
        for j = 1 : num
            sub_line = fgetl(fr); 
            space_idx = isspace(sub_line); 
            data(j) = str2num(sub_line(1:find(space_idx==1, 1)-1)); 
        end
        dlmwrite([root_file 'Ristretto_layer_' num2str(layer_cnt), '_w.csv'], data);
        fprintf('Finish Processing Layer %d\n', layer_cnt);
    end
    
end