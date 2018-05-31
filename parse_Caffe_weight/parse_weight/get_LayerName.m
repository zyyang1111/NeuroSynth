function [layer_name] = get_LayerName(line)

space_idx = isspace(line); 
flag = 0;
for i = 1 : length(line)
    if space_idx(i) && ~space_idx(i+1)
        s = i+1; 
        flag = 1; 
    end
    if ~space_idx(i) && space_idx(i+1)
        if flag
            e = i;
            break
        end
    end
end

layer_name = line(s : e);