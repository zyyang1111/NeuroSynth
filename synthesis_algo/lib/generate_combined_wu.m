function [w_u] = generate_combined_wu(root_file, w2suffix)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zhiyuan Yang    June 29 2017
% Generate the intersected w_u 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[fig_size, file_array] = determine_NetInfo('AlexNet_caffe');
[fig_size, file_array] = read_NetInfo(root_file);
     
 idx = 1; 
 file = [root_file file_array{idx} '_w2' w2suffix '.csv'];
 w = csvread(file);
 w = w(:);
 
% Step 1: all the multiplications are separated; find the unique weights
% and the number of such multiplications 
w_u = unique(w); 
for i = 1 : size(w_u, 1)
    w_u(i,2) = length(find(w == w_u(i,1))) * fig_size(idx)^2;
end

w_u_re = w_u; 

for idx = 2 : length(fig_size)
    fprintf('Layer %d\n', idx);
    size_wu = size(w_u, 1);
    file = [root_file file_array{idx} '_w2' w2suffix '.csv'];
    w = csvread(file);
    w = w(:);
    w_u_t = unique(w); 
    ct = 0; 
    for i = 1 : size(w_u_t, 1)
        jdx = find(w_u(:,1) == w_u_t(i));
        if isempty(jdx)
            ct = ct + 1; 
            w_u(size_wu+ct, 1) = w_u_t(i); 
            w_u(size_wu+ct, 2) = length(find(w == w_u_t(i,1))) * fig_size(idx)^2;
        else
            w_u(jdx, 2) = max(w_u(jdx, 2), length(find(w == w_u_t(i,1))) * fig_size(idx)^2);
        end
    end
end

fprintf('Original size of w_u is %d\n', sum(w_u_re(:,2)));
fprintf('New size of w_u is %d\n', sum(w_u(:,2)));
fprintf('Ratio between the two sizes is %.4f\n', sum(w_u(:,2))/sum(w_u_re(:,2)));
