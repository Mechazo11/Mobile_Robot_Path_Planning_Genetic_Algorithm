function rand_ls = random_generator(num_to_gen, min_lim, max_lim) 
% This function by default will generate random numbers between min_lim and
% max_lim
    %rng(0.005,'twister');
    rng shuffle
    a = min_lim;
    b = max_lim;
    rand_ls = (b-a).*rand(num_to_gen,1) + a; % This is column vector
    rand_ls = transpose(rand_ls); % Swap up to row vector
end