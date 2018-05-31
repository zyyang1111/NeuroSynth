function [time, cycle_tot, cycle_compute, mem_in, mem_ke, mem_out, cycle_mem] = cycle_estimation_funct(Tn0, Tr0, M, N, R, C, K, S, Tm, Naddin, maxK, cache_line, mem_lat, max_mult_var, clk, conv)


if conv
    Tn = min(Tn0, N);
    Tr = min(Tr0, R); 

    % estimate:
    cycle_compute = ceil(M/Tm) * (fix(N./Tn) .* (R*C * (ceil(Tn * K^2 / Naddin))) + R*C*(ceil((N - fix(N./Tn).*Tn) * K^2 / Naddin)));
%     fprintf('Single simulation cycle : %d\n', (ceil(Tn * K^2 / Naddin)));

    mem_in = ceil(M/Tm) * N .* ( ((R-1)*S+K)^2 + fix((R-1)./Tr)*max(0,K-S)*((R-1)*S+K)*2 + (max(0,K-S))^2*(fix((R-1)./Tr)).^2 ); 
    mem_ke = (ceil(R./Tr)).^2 * M * N * K^2; 
    mem_out = M * R * C; 
    cycle_mem = ceil(((mem_in + mem_out) * 2 + mem_ke * max_mult_var / 8) / cache_line) * mem_lat; 

% cycle_stall = ceil(N./Tn) .* ceil(M/Tm) .* R .* C .* (ceil(Tn * K^2 / Naddin) - ceil(Tn*K^2/num_mult)); 

    cycle_tot = cycle_mem + cycle_compute;
    time = cycle_tot' * clk * 1e-6; % [ms]
else
    Tn = Tn0 * maxK^2; 
    Tn = min(Tn, N); 
    
    cycle_compute = ceil(M/Tm) * (fix(N./Tn) .* ceil(Tn / Naddin) + (ceil((N - fix(N./Tn).*Tn) / Naddin)));
    mem_in = ceil(M/Tm) * N; 
    mem_ke = ceil(M/Tm) * N * Tm; 
    mem_out = M; 
    cycle_mem = ceil(((mem_in + mem_out) * 2 + mem_ke * max_mult_var / 8) / cache_line) * mem_lat; 
    cycle_stall = ceil(Tn/Naddin) .* 0.5 .* ceil(M/Tm) .* fix(N./Tn);
    cycle_mem = cycle_mem - cycle_stall;
    
    cycle_tot = cycle_mem + cycle_compute;
    time = cycle_tot' * clk * 1e-6; % [ms]
end