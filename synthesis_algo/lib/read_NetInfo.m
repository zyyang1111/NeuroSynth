function [fig_size, file_array] = read_NetInfo(root_file)

% Zhiyuan Yang  June 8 2017 
% This function read the names of the layer and the size of the output
% feature map given the root_file

fr = fopen([root_file 'size.txt']);

idx = 0; 
while ~feof(fr)
    line = fgetl(fr); 
    if line(1) ~= '('
        idx = idx + 1; 
        bs = find(line == '/'); 
        if ~isempty(bs)
            line(bs) = '_';
        end
        file_array{idx} = line; 
    else
        num_dim = length(find(line == ','))+1; 
        if num_dim ~= 4
            fig_size(idx) = 1; 
        else
            b_s = find(line == ',') + 2;
            b_e = find(line == ')') - 1; 
            dim = str2double(line(b_s(end):b_e)); 
            fig_size(idx) = dim;
        end
    end
end

fig_size = fig_size'; 
file_array = file_array';

