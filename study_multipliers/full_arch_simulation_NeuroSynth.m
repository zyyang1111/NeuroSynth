% This script is built based on full_arch_simulation_V3.m

close all
clear
clc

addpath '../synthesis_algo/';
addpath '../synthesis_algo/num_convert/';
addpath './lib/';

root_root_file = '../../AlexNet_caffe/'; 
Q = 7;
%%%%%%%%%%%%%%%% Neural Network Info %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fig_size stores the size of the output figure
[fig_size, file_array] = read_NetInfo(root_root_file);
% % % AlexNet
adjust_size = 0;
if exist([root_root_file 'size_info.mat'])
    load([root_root_file 'size_info.mat']);
else
    [strid_a, M_a, N_a, K_a, layer_type] = read_otherNetInfo_new(root_root_file, file_array, adjust_size);
end


% strid_a = [4 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1];
% M_a = [96 ; 256 ; 384 ; 384 ; 256 ; 4096 ; 4096 ; 1000]; % number of output features
% N_a = [3 ; 96 ; 256 ; 384 ; 384 ; 9216 ; 4096 ; 4096]; % number of input features
% K_a = [11 ; 5 ; 3 ; 3 ; 3 ; 1 ; 1 ; 1];
num_layer = length(fig_size);

multQ_scale = 1;
addQ_scale = 6;


for l = 1 : num_layer
% for l = 4
    fprintf('============== LAYER %d =============\n', l);
% l = 4; % the layer in the NN to be studied 
R = fig_size(l); C = R;
S = strid_a(l); K = K_a(l);
M = M_a(l); N = N_a(l);

%%%%%%%%%%%%%%%%% Multipliers Info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the multipliers
area_limit = 0.01 * 1e8;
str1 = '001'; 
str2 = '001';
w2suffix = '';
w2suffix_2 = '';
suffix = '_comb';
sample_layer = 0;
root_file = [root_root_file 'separate_results/P' str1 '_A' str2 suffix w2suffix '/total_results/']; 
load([root_file 'multipliers.mat']);
% % modify multipliers
multipliers = modify_multipliers(multipliers, area_limit); 
% multipliers = modify_multipliers_round(multipliers, root_root_file, w2suffix_2);
%%%%%%%%%%%%%%%%%% Weight Info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = csvread([root_root_file file_array{l} '_w2' w2suffix_2 '.csv']);
w = round(w * 2^Q);
scale = N / (size(w,1) / K^2); 
if adjust_size
    w0 = w; 
    while scale > 1
        w0 = [w0 ; w];
        scale = scale - 1; 
    end
    w = w0;
end
% mapping each weight to the unique multiplier
w_map = map_w2mult(w, multipliers); 

%Tmm_a = 80 : -2 : 8;
%Trr_a = 1 : 1 : 100; 
%Naddin = 2:2:100; 

Tmm_a = 74; 
Trr_a = 8; 
Naddin = 16;
% Trr_a = 8 : 4 : 100;
% Trr_a = [16]; 
%sim_case = zeros(length(Trr_a), 4);
area_const = 5 * 1e6; 
cnt = 0;
for iNadd = 1 : length(Naddin)
    num_add_input = Naddin(iNadd);
for iTmmx = 1 : length(Tmm_a)
    Tmm = Tmm_a(iTmmx);
for iTrrx = 1 : length(Trr_a)
    Trr = Trr_a(iTrrx);
    Tcc = Trr; 
    flag = 0; tot_area = 0;
    Tnn = 0; 
    while tot_area <= area_const
        Tnn = Tnn + 2; 

%%%%%%%%%%%%%%%%%%% Control Info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tr = 1024; Tc = Tr; % output tile for L2 Cache 
% Trr = 170; Tcc = Trr; % output tile for L1 Cache 
Tm = 1024;
Tn = 1024; 
% Tmm = 45                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                ; 
% Tnn = 1; 
% num_add_input = 50; 
Outbuff_deep = Trr * Tcc; 
% cache_line = 128;
% define the queue-size
mult_queue_size = multipliers(:,3) * multQ_scale;
add_queue_size = num_add_input * addQ_scale;
%%%%%%%%%%%%%%%%%%% Hardware Info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------- Energy ---------------------------
% adder energy [nJ]
add_delay_0 = 1702; % ps
add_power_0 = 0.03755; % mW
add_energy_0 = add_delay_0 * add_power_0 * 1e-15 * 1e9; % nJ
% -------------- Memory Latency -------------------
% memory latency 
lat_main = 0; 
lat_L2 = 100;
mem_band = 250; %B/ns
% -------------- Area -----------------------------
[mult_lib, minW, maxW, len_mult_lib] = generate_mullib(); 
add_area_0 = 661; % um^2 
% mem_area_0 = 10.65; % um^2, 90nm tech
mem_area_0 = 54.6; % um^2, according to DianNao
% mem_area_0 = 208.9; % um^2


%%%%%%%%%%% --------------- estimate size -------------- %%%%%%%%%%
max_L2_out_size = Tr * Tc * Tm * 2; 
max_L2_in_size = max(((Tr-1)*strid_a+K_a) .* ((Tc-1)*strid_a+K_a) * Tn * 2); 
max_L2_ke_size = max(ceil(M_a .* N_a .* K_a.^2 * max(multipliers(:,2)) / 8));
max_L2_cache_size = (max_L2_out_size + max_L2_in_size + max_L2_ke_size) / (1024^2); % in MB
% fprintf('Required L2 cache size = %.3f (MB)\n', max_L2_cache_size);
% max_L1_out_size = Trr * Tcc * Tmm * 2; 
% max_L1_in_size = Trr * Tcc * Tnn * max(K_a)^2 * 2; 
% max_L1_ke_size = ceil(Tmm * Tnn * max(K_a)^2 * max(multipliers(:,2)) / 8); 
% max_L1_cache_size = (max_L1_out_size + max_L1_in_size + max_L1_ke_size) / (1024); % in KB 


max_L1_out_size = Trr * Tcc * Tmm * 2; 
max_L1_in_size = max(((Trr-1)*strid_a+K_a) .* ((Tcc-1)*strid_a+K_a) * Tnn * 2); 
max_L1_ke_size = ceil(Tmm * Tnn * max(K_a)^2 * max(multipliers(:,2)) / 8); 
max_L1_cache_size = (max_L1_out_size + max_L1_in_size + max_L1_ke_size) / (1024); % in KB 
% fprintf('Required L1 cache size = %.3f (KB)\n', max_L1_cache_size);
%%%%%%%%%%% --------------- estimate area -------------- %%%%%%%%%%
mult_area = sum((mult_lib(multipliers(:,1)-minW+1, 3) + len_mult_lib(multipliers(:,2)+1, 2)) .* multipliers(:,3));
add_area = (num_add_input - 1) * Tmm * add_area_0; 
mem_area = (max_L1_cache_size * 1024 + sum(mult_queue_size) * 2 + add_queue_size * Tmm * 2) * mem_area_0; 
tot_area = mult_area + add_area + mem_area;
% fprintf('Required Chip Area = %.3f (mm^2)\n', tot_area / 1e6);
if tot_area < area_const
    flag = 1; 
end
    end
    if ~flag
        break;
    end
    Tnn = Tnn - 2;
     
    %%%%%%%%%%% --------------- estimate size -------------- %%%%%%%%%%
    max_L2_out_size = Tr * Tc * Tm * 2; 
    max_L2_in_size = max(((Tr-1)*strid_a+K_a) .* ((Tc-1)*strid_a+K_a) * Tn * 2); 
    max_L2_ke_size = max(ceil(M_a .* N_a .* K_a.^2 * max(multipliers(:,2)) / 8));
    max_L2_cache_size = (max_L2_out_size + max_L2_in_size + max_L2_ke_size) / (1024^2); % in MB
    fprintf('Required L2 cache size = %.3f (MB)\n', max_L2_cache_size);
%     max_L1_out_size = Trr * Tcc * Tmm * 2; 
%     max_L1_in_size = Trr * Tcc * Tnn * max(K_a)^2 * 2; 
%     max_L1_ke_size = ceil(Tmm * Tnn * max(K_a)^2 * max(multipliers(:,2)) / 8); 
%     max_L1_cache_size = (max_L1_out_size + max_L1_in_size + max_L1_ke_size) / (1024); % in KB 
    max_L1_out_size = Trr * Tcc * Tmm * 2; 
    max_L1_in_size = max(((Trr-1)*strid_a+K_a) .* ((Tcc-1)*strid_a+K_a) * Tnn * 2); 
    max_L1_ke_size = ceil(Tmm * Tnn * max(K_a)^2 * max(multipliers(:,2)) / 8); 
    max_L1_cache_size = (max_L1_out_size + max_L1_in_size + max_L1_ke_size) / (1024); % in KB 
    
    fprintf('Required L1 cache size = %.3f (KB)\n', max_L1_cache_size);
    %%%%%%%%%%% --------------- estimate area -------------- %%%%%%%%%%
    mult_area = sum((mult_lib(multipliers(:,1)-minW+1, 3) + len_mult_lib(multipliers(:,2)+1, 2)) .* multipliers(:,3));
    add_area = (num_add_input - 1) * Tmm * add_area_0; 
    mem_area = (max_L1_cache_size * 1024 + sum(mult_queue_size) * 2 + add_queue_size * Tmm * 2) * mem_area_0; 
    tot_area = mult_area + add_area + mem_area;
    fprintf('Required Chip Area = %.3f (mm^2)\n', tot_area / 1e6);
    fprintf('Tnn = %d\n', Tnn);
    cnt = cnt + 1; 
   



%%%%%%%%%%% --------------- estimate clk -------------- %%%%%%%%%%
mult_clk = max(max(mult_lib(multipliers(:,1)-minW+1, 2), len_mult_lib(multipliers(:,2)+1, 1)));
add_clk = add_delay_0 * ceil(log(num_add_input) / log(2)); 
% % % % determine using the slowest module
clk = max(mult_clk , add_clk) / 1e3; % ns
clk_n = 1;
% % % % determine according to the multiplier
% if add_clk <= mult_clk
%     clk = mult_clk / 1e3; 
% else
%     clk_n = 1; 
%     while add_clk / clk_n > mult_clk
%         clk_n = clk_n + 1; 
%     end
%     clk_n = clk_n - 1;
%     clk = add_clk / clk_n / 1e3;
% end

if l == 1
    sim_case(cnt, 1:8) = [Tmm, num_add_input, Trr, Tnn, multQ_scale, addQ_scale, tot_area/1e6, clk]; 
end

cache_line = fix(clk * mem_band);

% adjust for FC layers 
if ~layer_type(l) % for AlexNet
    Tm = Tm * Tr * Tc; % Input Size for the L2 cache
    Tn = max(((Tr-1)*S+K_a) .* ((Tc-1)*S+K_a) * Tn); 
%     Tnn = Trr * Tcc * Tnn * max(K_a)^2; 
    Tnn = Tnn * 50;
    Tr = 1; Tc = 1; 
    Trr = 1; Tcc = 1;
end

% % build multiplier queue 
mult_queue = zeros(size(multipliers, 1), 1); % stores the number of unprocessed inputs 
mult_out_ind = zeros(size(multipliers, 1), 1); % stores the index of the output 
mult_out_num = zeros(size(multipliers, 1), 1); % stores the number of operations for that output
% % build adder queue
add_queue = zeros(Tmm, 1); % stores the number of unprocessed addition oprand
add_out_ind = zeros(Tmm, 1); % stores the index of the output
add_out_num = zeros(Tmm, 1); % stores the number of operands for the output
% % build output queue
out_done = [];
out_done_num = size(w, 1); % number of adding oprand to finish a output
out_ind_base_lib = [];

%%%%%%%%%%%%%%%%%%%% Main Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L1_in = 0;
L1_ke = 0;
L1_ke_2 = 0;
L1_out = 0;
L1_stall = 0;

cycle_num = 0;
tot_energy = 0;
stall_cycle = 0;
add_energy = 0;
mult_energy = 0;

main_mem_cycle = 0; 
main_mem_energy = 0;
L2_mem_cycle = 0;
L2_mem_energy = 0; 
% memory access statistics
load_from_main = 0; 
load_from_L2 = 0;
write_to_L2 = 0; 
write_to_main = 0;
% cache tags
L2_out_tag = -1; 
L1_out_tag = -1;


max_real_L1_ke_size = 0;
max_real_L1_in_size = 0;
max_real_L1_out_size = 0;
% if it is the convolutional layer
L2_ke_size = ceil(M * N * K^2 * max(multipliers(:,2)) / 8); % in bytes 
% fprintf('Load L2 kernel Size is %d B\n', L2_ke_size);
load_from_main = load_from_main + L2_ke_size; 
for row = 1 : Tr : R
    tr = min(Tr, R-row+1); 
    for col = 1 : Tc : C
        tc = min(Tc, C-col+1);
        for to = 1 : Tm : M
            tm = min(Tm, M-to+1);
            
            % load output from main space 
            [L2_out_size, L2_out_dim] = L2_declare_output(row, col, to, tr, tc, tm); 
            L2_out_tag_n = (L2_out_dim(1,1)*100+L2_out_dim(2,1)) * 1000 + L2_out_dim(3,1); 
            if L2_out_tag_n ~= L2_out_tag
                if L2_out_tag > 0
                    % write output back to the main memory 
                    % fprintf('\tStore Out to Main (%d B) for (%d, %d, %d) to (%d, %d, %d)\n', ...
                        % L2_to_main_out_size, L2_to_main_out_dim(1,1), L2_to_main_out_dim(2,1), L2_to_main_out_dim(3,1), ...
                        % L2_to_main_out_dim(1,2), L2_to_main_out_dim(2,2), L2_to_main_out_dim(3,2)); 
                    write_to_main = write_to_main + L2_to_main_out_size;
                    % write-to-main latency and energy 
                    cycle_num = cycle_num + lat_main; main_mem_cycle = main_mem_cycle + lat_main;
                    tot_energy = tot_energy + getMainEnergy(L2_to_main_out_size); main_mem_energy = main_mem_energy + getMainEnergy(L2_to_main_out_size);
                end
                L2_out_tag = L2_out_tag_n;
                % update output load statistics
                % fprintf('Load L2 Output (%d B) for (%d, %d, %d) to (%d, %d, %d)\n', ...
                    % L2_out_size, L2_out_dim(1,1), L2_out_dim(2,1), L2_out_dim(3,1), ...
                    % L2_out_dim(1,2), L2_out_dim(2,2), L2_out_dim(3,2));
                load_from_main = load_from_main + L2_out_size; 
                % load-from-main out latency and energy
                tot_energy = tot_energy + getMainEnergy(L2_out_size); 
                main_mem_energy = main_mem_energy + getMainEnergy(L2_out_size);
            end

            for ti = 1 : Tn : N
                tn = min(Tn, N-ti+1);
                % L2-cache load 
                [L2_in_size, L2_in_tag] = L2_load_input(row, col, ti, tr, tc, tn, S, K); 
                % fprintf('Load L2 Input (%d B) for (%d, %d, %d) to (%d, %d, %d)\n', ...
                    % L2_in_size, L2_in_tag(1,1), L2_in_tag(2,1), L2_in_tag(3,1), ...
                    % L2_in_tag(1,2), L2_in_tag(2,2), L2_in_tag(3,2));
                
                cur_L2_cache_size = L2_ke_size + L2_out_size + L2_in_size; 
                % convert cur_L2_cache_size to MB
                cur_L2_cache_size = cur_L2_cache_size / 1024 / 1024; 
                % fprintf('Total used L2 cache size is %.3f (MB)\n', cur_L2_cache_size);
                
                % load-from-main in latency and energy
                cycle_num = cycle_num + lat_main; main_mem_cycle = main_mem_cycle + lat_main;
                tot_energy = tot_energy + getMainEnergy(L2_in_size); main_mem_energy = main_mem_energy + getMainEnergy(L2_in_size);
                
                load_from_main = load_from_main + L2_in_size;
                
                % go to on-chip processing 
                for rr = row : Trr : row+tr-1
                    trr = min(Trr, row+tr-rr);
                    for cc = col : Tcc : col+tc-1
                        tcc = min(Tcc, col+tc-cc);
                        for too = to : Tmm : to+tm-1
                            tmm = min(Tmm, to+tm-too); 
                            % load output from L2 to L1
                            L2_to_L1_out_size = trr * tcc * tmm * 2;
                            L1_out_tag_n = (rr*100+cc)*1000 + too;
                            if L1_out_tag_n ~= L1_out_tag
                                if L1_out_tag > 0
                                    % write output to L2-cache in bytes 
                                    % L1_to_L2_out_size is calculated before
%                                   L1_to_L2_out_size = trr * tcc * tmm * 2; 
                                    % fprintf('\tStore Out to L2 (%d B) for (%d, %d, %d) to (%d, %d, %d)\n', ...
                                        % L1_to_L2_out_size, L1_to_L2_out_dim(1,1), L1_to_L2_out_dim(2,1), L1_to_L2_out_dim(3,1), ...
                                        % L1_to_L2_out_dim(1,2), L1_to_L2_out_dim(2,2), L1_to_L2_out_dim(3,2)); 
                            
                                    write_to_L2 = write_to_L2 + L1_to_L2_out_size;
                                    % write-to-L2 out latency and energy
%                                     cycle_num = cycle_num + lat_L2; L2_mem_cycle = L2_mem_cycle + lat_L2;
                                    tot_energy = tot_energy + getL2Energy(L1_to_L2_out_size, 'w'); L2_mem_energy = L2_mem_energy + getL2Energy(L1_to_L2_out_size, 'w');
                                end
                                L1_out_tag = L1_out_tag_n; 
                                % update output load statistics
                                % fprintf('\tLoad L1 Output (%d B) for (%d, %d, %d) to (%d, %d, %d)\n', ...
                                    % L2_to_L1_out_size, rr, cc, too, rr+trr-1, cc+tcc-1, too+tmm-1); 
                                load_from_L2 = load_from_L2 + L2_to_L1_out_size; 
                                if L2_to_L1_out_size > max_real_L1_out_size
                                    max_real_L1_out_size = L2_to_L1_out_size;
                                end
                                % load-from-L2 out latency and energy
                                cycle_num = cycle_num + lat_L2 * ceil(L2_to_L1_out_size / cache_line); L2_mem_cycle = L2_mem_cycle + lat_L2 * ceil(L2_to_L1_out_size / cache_line);
                                tot_energy = tot_energy + getL2Energy(L2_to_L1_out_size, 'r'); 
                                L2_mem_energy = L2_mem_energy + getL2Energy(L2_to_L1_out_size, 'r');
                                L1_out = L1_out + L2_to_L1_out_size; 
                            end
                            
                            for tii = ti : Tnn : ti+tn-1
                                tnn = min(Tnn, ti+tn-tii);
                                % load input from L2 to L1
                                [L1_store_in_size, L1_load_in_size, L1_in_tag] = ...
                                    L1_load_input(rr, cc, tii, trr, tcc, tnn, S, K); 
                                if L1_store_in_size > max_real_L1_in_size
                                    max_real_L1_in_size = L1_store_in_size;
                                end
%                                 fprintf('\tLoad L1 Input (%d B) for (%d, %d, %d) to (%d, %d, %d) and stored in %d B space\n', ...
%                                     L1_load_in_size, L1_in_tag(1,1), L1_in_tag(2,1), L1_in_tag(3,1), ...
%                                     L1_in_tag(1,2), L1_in_tag(2,2), L1_in_tag(3,2), L1_store_in_size);
                                % load kernel from L2 to L1
                                L1_ke_size = ceil(tmm * tnn * K^2 * max(multipliers(:,2)) / 8); 
                                
                                % find the weight mapped to the
                                % single-weight multiplier
                                weight_prop = w((tii-1)*K*K+1 : (tii+tnn-1)*K*K, too:too+tmm-1);
                                wu = unique(weight_prop); 
                                L1_ke_size_2 = 0;
                                for iwu = 1 : length(wu)
                                    wu_idx = find(w_map(:,1) == wu(iwu)); 
                                    if multipliers(w_map(wu_idx, 2), 2) > 0
                                        L1_ke_size_2 = L1_ke_size_2 + length(find(weight_prop == wu(iwu))); 
                                    end
                                end
                                L1_ke_size_2 = ceil(L1_ke_size_2 * max(multipliers(:,2)) / 8); 
                                L1_ke_size = L1_ke_size_2; 
                                if L1_ke_size > max_real_L1_ke_size
                                    max_real_L1_ke_size = L1_ke_size; 
                                end
%                                 fprintf('\tLoad L1 kernel Size is %d B\n', L1_ke_size);
                                
                                % load-from-L2 in and kernel latency and
                                % energy
%                                 cycle_num = cycle_num + lat_L2 * max(0, ceil((L1_load_in_size + L1_ke_size)/ cache_line)-stall_cycle); L2_mem_cycle = L2_mem_cycle + lat_L2 * max(0, ceil((L1_load_in_size + L1_ke_size)/ cache_line)-stall_cycle); 
                                cycle_num = cycle_num + max(0, lat_L2 * (ceil((L1_load_in_size + L1_ke_size)/ cache_line)) - stall_cycle); 
                                L2_mem_cycle = L2_mem_cycle + max(0, lat_L2 * (ceil((L1_load_in_size + L1_ke_size)/ cache_line)) - stall_cycle); 
                                
                                L1_stall = L1_stall + stall_cycle; 
                                
                                tot_energy = tot_energy + getL2Energy(L1_load_in_size + L1_ke_size, 'r'); L2_mem_energy = L2_mem_energy + getL2Energy(L1_load_in_size + L1_ke_size, 'r');
                                
                                cur_L1_cache_size = L1_store_in_size + L2_to_L1_out_size + L1_ke_size; 
                                % convert cur_L1_cache_size to KB
                                cur_L1_cache_size = cur_L1_cache_size / 1024; 
                                % fprintf('\tTotal used L1 cache size is %.3f (KB)\n', cur_L1_cache_size);
                        
                                load_from_L2 = load_from_L2 + L1_load_in_size + L1_ke_size; 
                                
                                L1_ke = L1_ke + L1_ke_size; 
                                L1_in = L1_in + L1_load_in_size; 
                                
                                % go to process compute unites
                                for itr = 1 : trr 
                                    for itc = 1 : tcc 
                                        % develop the out_done array
                                        out_ind_base = (((rr+itr-1)*100 + (cc+itc-1))*10 + fix(too/Tmm))*100;
                                        % out_ind_base is used to generate
                                        % a unique ID for each output
                                        if isempty(find(out_ind_base_lib == out_ind_base, 1))
                                            out_done_a = zeros(tmm, 2); 
                                            for ii = 1 : size(out_done_a, 1)
                                                out_done_a(ii, 1) = ii + out_ind_base;
                                            end
                                            % yzy: we can also preallocate
                                            % the space for the following
                                            % two arrays
                                            out_done = [out_done ; out_done_a];
                                            out_ind_base_lib = [out_ind_base_lib ; out_ind_base];
                                        end
                                        
                                        % process one cache-line 
                                        [mult_use] = gen_multuse(too, tii, tmm, tnn, K, w);
                                        
                                        
                                        cycle_0 = cycle_num; 
                                        [mult_queue, mult_out_ind, mult_out_num, add_queue, add_out_ind, add_out_num, out_done, cycle_num, tot_energy, stall_cycle, mult_energy, add_energy] = ...
                                            one_cycle_pipeline_4(multipliers, mult_use, mult_queue, mult_out_ind, mult_out_num, add_queue, add_out_ind, add_out_num, ...
                                            out_ind_base, tmm, num_add_input, out_done, out_done_num, cycle_num, mult_queue_size, add_queue_size, add_energy_0, ...
                                            tot_energy, w_map, clk_n, add_energy, mult_energy);
                                        
%                                         fprintf('\tStall Cycles is %d\n', stall_cycle);
%                                         fprintf('Total Cycles=%d\tIncreased Cycles is %d\n', cycle_num, cycle_num - cycle_0);
                                    end % itc
                                end % itr

                            end % tii
                            if length(out_done) >= L2_to_L1_out_size / 2
                                % if the output buffer is full
                                % stall until the output buffer is empty
                                [add_queue, add_out_ind, add_out_num, out_done, tot_energy] = empty_outbuffer(add_queue, add_out_ind, add_out_num, num_add_input, out_done, out_done_num, add_queue_size, add_energy_0, tot_energy);
                            end
                            % record the current L1 output information
                            L1_to_L2_out_size = trr * tcc * tmm * 2; 
                            L1_to_L2_out_dim = [rr, rr+trr-1; cc, cc+tcc-1; too, too+tmm-1];
                        end % too
                    end % cc
                end % rr
            end % ti
            % record the current L2 output information
            L2_to_main_out_size = L2_out_size; 
            L2_to_main_out_dim = L2_out_dim;
            
        end % to
    end % col
end % row

extra_cycle_num = 0;
while max(mult_queue) ~= 0 || max(add_queue) ~= 0
    cycle_num = cycle_num + 1;
    extra_cycle_num = extra_cycle_num + 1;
%     [mult_queue, mult_out_ind, mult_out_num, add_queue, add_out_ind, add_out_num, tot_energy] = ...
%         process_mult_new3(multipliers, mult_queue, mult_out_ind, mult_out_num, add_queue, add_out_ind, add_out_num, Tm, tot_energy);
    [mult_queue, mult_out_ind, mult_out_num, add_queue, add_out_ind, add_out_num, tot_energy] = ...
        process_mult_new2(multipliers, mult_queue, mult_out_ind, mult_out_num, add_queue, add_out_ind, add_out_num, tm, tot_energy, out_ind_base);
    if max(mult_out_ind(:,1))+1 < size(mult_out_ind, 2)
        mult_out_ind(:, max(mult_out_ind(:,1))+2:end) = [];
        mult_out_num(:, max(mult_out_ind(:,1))+2:end) = [];
    end
    [add_queue, add_out_ind, add_out_num, out_done, tot_energy] = process_add(add_queue, add_out_ind, add_out_num, num_add_input, out_done, out_done_num, add_queue_size, add_energy_0, tot_energy);
    if max(add_out_ind(:,1))+1 < size(add_out_ind, 2)
        add_out_ind(:, max(add_out_ind(:,1))+2:end) = [];
        add_out_num(:, max(add_out_ind(:,1))+2:end) = [];
    end
    out_done(out_done(:,2) == out_done_num, :) = [];
end


exe_time = cycle_num * clk; 
fprintf('Total Execution Time is %.3f [ms]\n', exe_time * 1e-6);
sim_case(cnt, 8+l) = exe_time * 1e-6; 
fprintf('Total Cycle %.0f\n', cycle_num);
fprintf('Compute Cycle %.0f\n', cycle_num - L2_mem_cycle);
fprintf('Mem Cycle %.0f\n', L2_mem_cycle);
fprintf('Total Energy: %.5f\n', tot_energy);
fprintf('[mult, add, mem] energy: %.5f, %.5f, %.5f\n', mult_energy, add_energy, tot_energy-mult_energy-add_energy);
 

results(:,l) = [exe_time * 1e-6; cycle_num; cycle_num - L2_mem_cycle; L2_mem_cycle; tot_energy; mult_energy; add_energy; tot_energy-mult_energy-add_energy;L1_in;L1_out;L1_ke];
end
end
end
end
% [L1_in, L1_out, L1_ke]
% 
% [max_real_L1_in_size, max_real_L1_out_size, max_real_L1_ke_size]
% 
% [max_L1_in_size, max_L1_out_size, max_L1_ke_size]