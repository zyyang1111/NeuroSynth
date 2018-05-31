function [strid_a, M_a, N_a, K_a] = read_otherNetInfo(root_file, file_array)

% Zhiyuan Yang  Sep 18 2017 
% This function reads the size of the K size
K_a = [];
strid_a = []; 
stride = []; 

fr = fopen([root_file 'deploy.prototxt']);
flag = 0;
while ~feof(fr)
    line = fgetl(fr); 
    
    type = strfind(line, 'type'); 
    if ~isempty(type)
        
        tp = strfind(line, 'Convolution');
        if ~isempty(tp)
            strid_a = [strid_a ; stride]; 
            flag = 1; 
            stride = 1;
        else
            flag = 0; 
            tp2 = strfind(line, 'InnerProduct');
            if ~isempty(tp2)
                strid_a = [strid_a ; stride]; 
                K = 1; 
                K_a = [K_a ; K];
                stride = 1; 
            end
        end
    end
    
    if flag
        kernel = strfind(line, 'kernel_size'); 
        if ~isempty(kernel)
            com = find(line == ':');
            K = str2num(line(com+1 : end));
            K_a = [K_a ; K];
        end
        
        strd = strfind(line, 'stride'); 
        if ~isempty(strd)
            com = find(line == ':');
            stride = str2num(line(com+1 : end));
        end
    end
end
strid_a = [strid_a ; stride]; 

% define M_a and N_a
M_a = []; 
N_a = []; 
for idx = 1 : length(file_array)
    K = K_a(idx);
    w = csvread([root_file file_array{idx} '_w2.csv']);
    N = size(w, 1) / K^2; 
    if idx > 1
        alpha = M / N; 
    else
        alpha = 1;
    end
    w0 = w; 
    while alpha > 1
        w0 = [w0 ; w]; 
        alpha = alpha - 1;
    end
    M = size(w0, 2);
    N = size(w0, 1) / K^2;
    
    M_a = [M_a ; M];
    N_a = [N_a ; N]; 
end

save([root_file 'size_info.mat'], 'strid_a', 'M_a', 'N_a', 'K_a'); 



