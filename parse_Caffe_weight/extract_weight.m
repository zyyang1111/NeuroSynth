% Zhiyuan Yang : May 20 2017 
% this script parses and read the weight from the 

close all 
clear
clc

addpath './parse_weight/';

root_file = '../VGG_16/'; 
fr = fopen([root_file 'new_res.txt'], 'r');

quantitize = 1;
quant_suffix = '';
suffix = ''; % '_round' or ''

if quantitize
    cluster = csvread([root_file 'weight_quant' quant_suffix '.csv']); 
    minW = cluster(1);
    num_cluster = cluster(2);
    LARGE = cluster(3); 
    cluster(1:3) = []; 
    suffix = ['_round' quant_suffix];
end

Q = 7;

while ~feof(fr)
    line = fgetl(fr); 
    
%     if isempty(line)
%         continue;
%     end
    
    if line(1) == '='
        layer_name = get_LayerName(line); 
        if ~isempty(find(layer_name == '/', 1))
            layer_name(layer_name == '/') = '_';
        end
        fprintf('======= %s ======\n', layer_name);
    end
    
    if line(1) == '['
        Ldim = get_Ldim(line);
        n_con = 1; 
        for i = 2 : length(Ldim)
            n_con = n_con * Ldim(i); 
        end
        weight = zeros(n_con, Ldim(1));
        weight_2 = zeros(n_con, Ldim(1));
        fprintf('['); 
        for i = 1 : length(Ldim)-1
            fprintf('%d, ', Ldim(i));
        end
        fprintf('%d]\n', Ldim(end));
    end
    
   if line(1) == 'C'
       nc = get_ConIdx(line); 
       if ~mod(nc, 50)
           fprintf('Connection %d\n', nc);
       end
       for i = 1 : n_con
           line = fgetl(fr); 
           dot = find(isspace(line) == 1); 
           if length(dot) > 1
               fprintf('More than 1 space!\n');
           end
           weight(i, nc) = str2double(line(1:dot-1));
           
           if quantitize
               % round the weight
               intweight = str2double(line(dot+1:end)); 
               if intweight >= minW && intweight < minW + num_cluster
                   if cluster(intweight - minW + 1) < LARGE
                       intweight = cluster(intweight - minW + 1); 
                   end
               end
               weight_2(i, nc) = intweight / 2^Q; 
           else
               weight_2(i,nc) = str2double(line(dot+1:end)) / 2^Q;
           end

           dum = fgetl(fr); % for GoogleNet
       end
        
       if nc == Ldim(1)
           % dlmwrite([root_file layer_name, '_w1.csv'], weight);
           dlmwrite([root_file layer_name, '_w2' suffix '.csv'], weight_2, 'precision', 9); 
       end 
   end
    
end

