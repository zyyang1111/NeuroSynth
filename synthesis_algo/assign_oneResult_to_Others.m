close all 
clear
clc

root_file = '../../GenderNet/'; 
% area = [0.01, 0.02, 0.04, 0.08, 0.1, 0.2, 0.4, 0.6];
area = 0.01;
file_str = '001';
w2suffix = '_round_2';
suffix = ['_comb' w2suffix];
% suffix = '';
sample_layer = 0;

power = zeros(length(area), 1);
time = zeros(length(area), 1);
clk = zeros(length(area), 1);
energy = zeros(length(area), 1); 

for i = 1 : length(area)
    str = num2str(area(i));
    fprintf('================ AREA (%s) ================\n', str);
    str(str == '.') = []; 
    
    [clk(i), time(i), energy(i), power(i)] = Multiplication_Mapping(root_file, file_str, str, suffix, w2suffix, sample_layer);
end

data = [time, energy, power];
