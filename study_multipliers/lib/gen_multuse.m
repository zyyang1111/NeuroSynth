function [mult_use] = gen_multuse(to, ti, tm, tn, K, w)

idx = 0;
mult_use = []; 
for m = to : to+tm-1
    idx = idx + 1; 
    set = w((ti-1)*K*K+1 : (ti+tn-1)*K*K, m);
    mult_use(idx).set = get_mult_use(set); 
    mult_use(idx).done = 0;
end