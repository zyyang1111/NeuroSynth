function [w] = gen_w(to, ti, tm, tn, K, w0)

idx = 0;
w = []; 
for m = to : to+tm-1
    idx = idx + 1; 
    set = w0((ti-1)*K*K+1 : (ti+tn-1)*K*K, m);
    w = [w ; set]; 
end