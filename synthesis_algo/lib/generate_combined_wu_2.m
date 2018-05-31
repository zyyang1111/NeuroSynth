function [w_u] = generate_combined_wu_2(root_file, w2suffix, Tm, Tn, Tr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zhiyuan Yang    Sep 18 2017
% Generate the intersected w_u 
% Different from "generate_combined_wu"
%   - consider the data access pattern 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[fig_size, file_array] = determine_NetInfo('AlexNet_caffe');
[fig_size, file_array] = read_NetInfo(root_file);
[K_a] = read_Ka(root_file);

cnt = 0; 

for idx = 1 : length(fig_size)
    fprintf('Layer %d\n', idx);
    file = [root_file file_array{idx} '_w2' w2suffix '.csv'];
    w = csvread(file);
    K = K_a(idx); 
    
    % # input feature maps
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
    % # output feature maps
    M = size(w, 2); 
    
    if K == 1
        Tn = Tn * Tr^2 * max(K_a)^2; 
    end
    
    
    for ti = 1 : Tn : N
        tn = min(Tn, N-ti+1);
        for to = 1 : Tm : M
            tm = min(Tm, M-to+1);
            w = gen_w(to, ti, tm, tn, K, w0); 
            w_u_t = unique(w); 
            if cnt == 0
                w_u = w_u_t; 
                for i = 1 : size(w_u, 1)
                    w_u(i,2) = length(find(w == w_u(i,1))) * fig_size(idx)^2;
                end
                w_u_re = w_u; 
                cnt = cnt + 1; 
            else
                size_wu = size(w_u, 1);
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
                cnt = cnt + 1; 
            end
        end
    end
end

fprintf('Original size of w_u is %d\n', sum(w_u_re(:,2)));
fprintf('New size of w_u is %d\n', sum(w_u(:,2)));
fprintf('Ratio between the two sizes is %.4f\n', sum(w_u(:,2))/sum(w_u_re(:,2)));
