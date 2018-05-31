function [K_a] = read_Ka(root_file)

% Zhiyuan Yang  Sep 18 2017 
% This function reads the size of the K size
K_a = [];
fr = fopen([root_file 'deploy.prototxt']);
flag = 0;
while ~feof(fr)
    line = fgetl(fr); 
    
    type = strfind(line, 'type'); 
    if ~isempty(type)
        tp = strfind(line, 'Convolution');
        if ~isempty(tp)
            flag = 1; 
        else
            flag = 0; 
            tp2 = strfind(line, 'InnerProduct');
            if ~isempty(tp2)
                K = 1; 
                K_a = [K_a ; K];
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
    end
end