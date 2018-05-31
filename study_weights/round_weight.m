close all
clear
clc

net_name = 'VGG_16'; 
net_name_suffix = '';
root_file = ['../../' net_name net_name_suffix '/']; 
addpath '../synthesis_algo/';
addpath '../synthesis_algo/num_convert/';
addpath './lib/';

w2suffix = '';

Q = 7; 
BIT_SIZE = 16;

[mult_lib, minW, maxW, len_mult_lib] = generate_mullib(); 
[fig_size, file_array] = read_NetInfo(root_file);
% file_array = {'conv1'; 'conv2'; 'ip1'; 'ip2'};
w_u = generate_combined_wu(root_file, w2suffix); % generate the weigth information of all layers
w_u(:,1) = w_u(:,1) * 2^Q; 
[~, sortID] = sort(w_u(:,1)); 
w_u = w_u(sortID, :);
fprintf('--------------------------------------------------\n');

% add the delay and area to "w_u"
for i = 1 : size(w_u, 1)
    w = w_u(i,1); 
    w_u(i, 3:4) = mult_lib(w-minW+1, 2:3); 
end

% use dynamic programming algorithm to find the optimal cluster
% the optimal cluster is the one which has both the smallest delay and area
% NOTE: we should determine the maximum size of the cluster 
cluster_size = 2;  
cluster_type = 0; % 1: don't consider wu_th; 0: consider wu_th
wu_th = 4; 
% % build hash table
% for i = 1 : size(w_u,1)
%     hash_tab(i).done = 0;
%     hash_tab(i).tot_D = 0;
%     hash_tab(i).tot_A = 0;
%     hash_tab(i).wc = []; 
% end
type = 1; % type = 1 --> round down 
          % type = 0 --> round to any value within the cluster
% w_u_n = w_u(abs(w_u(:,1)) <= wu_th, :); 
% w_u_n_idx = find(abs(w_u(:,1)) <= wu_th, 1) - 1; 
if cluster_type
    hash_tab = build_hash_tab(w_u);
    [tot_D, tot_A, wc, print_root, hash_tab] = determine_wc(w_u, cluster_size, 1, 1000, hash_tab, type);
else
    w_u_n1 = w_u(w_u(:,1) < -wu_th, :); 
    hash_tab = build_hash_tab(w_u_n1);
    [~, ~, wc1, ~, ~] = determine_wc(w_u_n1, cluster_size, 1, 1000, hash_tab, type); 
    w_u_n2 = w_u(w_u(:,1) > wu_th, :); 
    hash_tab = build_hash_tab(w_u_n2);
%     w_u_n2_idx = w_u(find(w_u(:,1) > wu_th, 1) - 1, 1); 
    w_u_n2_idx = find(w_u(:,1) > wu_th, 1) - 1;
    [tot_D, tot_A, wc2, print_root, hash_tab] = determine_wc(w_u_n2, cluster_size, 1, 1000, hash_tab, type); 
    wc2 = wc2 + w_u_n2_idx; 
    wc = [wc1' ; find(abs(w_u(:,1)) <= wu_th) ; wc2']; 
end
% wc = wc + w_u_n_idx; 
% wc = [find(w_u(:,1) < -wu_th); wc'; find(w_u(:,1) > wu_th)]; 
% wc = [wc' ; find(w_u(:,1) > wu_th)];
if type
    for i = 1 : length(wc)-1
        array = [];
        for j = wc(i) : size(w_u, 1)
            if j < wc(i+1)
                array = cat(2, array, j); 
            else
                break;
            end
        end
        cluster(i).w = wc(i);
        cluster(i).covA = array; 
    end
    i = length(wc);
    array = wc(i) : size(w_u, 1);
    cluster(i).w = wc(i);
    cluster(i).covA = array;
else
    for i = 1 : length(wc)
        cluster(i) = wc(length(wc)-i+1);
    end
end

% print to the file 
LARGE = 10000; 
minW = min(w_u(:,1)); 
maxW = max(w_u(:,1));
sizeW = maxW - minW + 1; 
data = LARGE * ones(sizeW, 1);
for i = 1 : length(cluster)
    tar = w_u(cluster(i).w, 1); 
    for j = 1 : length(cluster(i).covA)
        id = w_u(cluster(i).covA(j),1);
        data(id-minW+1) = tar; 
    end
end

% write to file 
if cluster_type
    filename = ['weight_quant_' num2str(cluster_size)];
else
    filename = ['weight_quant_' num2str(cluster_size) '_' num2str(wu_th)];
end
fw = fopen([root_file filename '_' net_name '.txt'], 'w');
fprintf(fw, '%d\n', min(w_u(:,1)));
fprintf(fw, '%d\n', length(data));
fprintf(fw, '%d\n', LARGE);
for i = 1 : length(data)
    fprintf(fw, '%d\n', data(i));
end
fclose(fw);

dlmwrite([root_file filename '.csv'], [minW; length(data); LARGE; data]);
            