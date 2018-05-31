function [tot_D, tot_A, wc, print_root, hash_tab] = determine_wc(w_u, cluster_size, root, print_root, hash_tab, type)

if hash_tab(root).done
    tot_D = hash_tab(root).tot_D;
    tot_A = hash_tab(root).tot_A;
    wc = hash_tab(root).wc;
else
    tot_D = inf;
    tot_A = inf; 
    wc = [];

    for i = 1 : cluster_size
        if i+1 > size(w_u, 1)
            if type
                D = calc_cost(w_u, i, 'D', cluster_size); 
                A = calc_cost(w_u, i, 'A', cluster_size);
                if D < tot_D && A < tot_A 
                    tot_D = D; tot_A = A; 
                    wc = root;
                end  
            else
                [D, ~] = calc_cost_new(w_u, i, 'D', cluster_size); 
                [A, sel, covA] = calc_cost_new(w_u, i, 'A', cluster_size);
                if D < tot_D && A < tot_A
                    tot_D = D; tot_A = A;
                    wc.w = sel+root-1;
                    wc.covA = covA+root-1;
                end
            end
            break
        else
            [D, A, wc_t, print_root, hash_tab] = determine_wc(w_u(i+1:end, :), cluster_size, root+i, print_root, hash_tab, type); 
            if type
                D = calc_cost(w_u, i, 'D', cluster_size) + D; 
                A = calc_cost(w_u, i, 'A', cluster_size) + A; 
                if D < tot_D && A < tot_A
                    tot_D = D; tot_A = A;
                    wc = [root, wc_t];
                end
            else
                D_t = calc_cost_new(w_u, i, 'D', cluster_size);  D = D + D_t; 
                [A_t, sel, covA] = calc_cost_new(w_u, i, 'A', cluster_size); A = A + A_t;
                if D < tot_D && A < tot_A
                    tot_D = D; tot_A = A;
                    wc = wc_t;
                    wc_len = length(wc);
                    wc(wc_len+1).w = sel+root-1;
                    wc(wc_len+1).covA = covA+root-1;
                end
            end
        end
    end
    hash_tab(root).done = 1; 
    hash_tab(root).tot_D = tot_D; 
    hash_tab(root).tot_A = tot_A; 
    hash_tab(root).wc = wc; 
end
if root < print_root
    fprintf('Finish Root %d\n', root);
    print_root = root;
end
        