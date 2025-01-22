function new_gen = create_new_gen(mmx,bit_count,point_ls)
    % HARDCODED, very specific to this problem
    % Generate a lookup table to convert integer positions into binary
    ax = point_ls(1,3); % min point index (will be 1)
    bx = point_ls(1,4); % max point index (i.e. 16,32 etc)
    loc_bunpack = [point_ls(1,1), point_ls(1,2)]; % Needed to preserve bits
    % in starting and ending positions
    
    % Make look up table
    % Initialize a matrix to hold binary array (max_point_index by bit_count) 
    look_up = load_lookup_table(bx,bit_count); % Shape bx by bit_count
    
    % Initialize working variables
    row_idx = size(mmx,1);
    col_idx = size(mmx,2);
    % Temporary matrix to hold children in N/2 by 2*n_feature shape
    % We will fill up new children path in this array
    cc_mat = zeros(row_idx, col_idx);
    
    % Arrays to hold children candidates in binary 
    ccx1 = [];
    ccx2 = [];
    
    % This part needs to be in a loop
    % Send one couple row, get two child path
    for ii = 1:row_idx
        %fprintf("\nCouple row selected %d",ii);
        in_mmx = mmx(ii,:); % One couple row
        [child1, child2] = path_genetic_cross(in_mmx, bit_count, look_up, loc_bunpack);
        % Combine two child paths into one set
        % Need to match N/2 * (2 * num_of_gene) shape in the end
        cons_mat = [child1, child2]; % Convert two row vector in 1
        % Push this children set into cc_mat matrix, gene by gene
        for k = 1: col_idx
            cc_mat(ii,k) = cons_mat(1,k);
        end
    end
    % Return all children path in appropriate shape
    new_gen = cc_mat; % shape --> (N/2 * (2* num_of_genes))
end

% ------------------------------------------------------------------------
function [child1,child2] = path_genetic_cross(in_mmx, bit_count, look_up, loc_bunpack)
    % Encode scheme
    % path-> bit_representation
    % 9 path gene -> @4bit, 72bit
    
    % Convert genes to binary, row vector
    in_mmx = DectoBin(in_mmx, bit_count, look_up); % Updated lookup
    bb_point = size(in_mmx,2)/2;
    
    % break cn*gene_num*bit_count*2 bit into two halves
    male_chromo = in_mmx(1,[1:bb_point]);
    female_chromo = in_mmx(1,[(bb_point+1):end]);
    
    % Preserve points
    first_lock = 0 + bit_count; %i.e @4 bit encoding, 1~4 position bit
    second_lock = length(male_chromo) - bit_count; 
    second_lock = second_lock + 1; % MATLAB does inclusive counting
    %i.e @4 bit encoding, 9 genes = 36bits. Thus, 36 -4 = 32 --> 32~36
    %position bit
    
    % Preserve starting and ending positions
    % HARDCODED, only works for 2D path finding problem
    % Use a function to populate p1 and p3 manually
    % This will effectively ensure the first and last gene never gets
    % changed
    [p1,p3] = preserve_pos(loc_bunpack, look_up);
    p2 = p1; % female start pos 
    p4 = p3; % female end pos
    
    % break out bits that can undergo cross over and mutation
    xx1 = (first_lock+1); % We have move one position to right
    xx2 = (second_lock-1); % We have to come one position to left
    male_x = male_chromo(1,[xx1:xx2]); % MATLAB does inclusive counting
    female_x = female_chromo(1,[xx1:xx2]);
    px_row = size(male_x,1);
    px_col = size(male_x,2);
    swap_possible_loc = px_col - 1;
    
    % Double cross over points (Original, randomized)
    rand_prob = random_generator(px_col, 0, 0.99);
    rand1 = randsample(rand_prob, 1);
    rand2 = randsample(rand_prob, 1);
    a1 = int16(swap_possible_loc * rand1);
    a2 = int16(swap_possible_loc * rand2);
    
    % Double cross over points, HARDCODED  
    %a1 = 5;
    %a2 = 23;
    
    % Rearrange two numbers in ascending order
    ax = sort([a1 a2]);
    cxx1 = ax(1,1);
    cxx2 = ax(1,2);
    
    % For single point crossover
    cxx_one = max(ax);
    
    % ----------------------- Crossover operators ----------------
    
    % Perform double point cross-over
    % cc1,cc2 are bits between xx1 and x22
    %cc1 = [male_x(1:cxx1) female_x(cxx1+1:cxx2) male_x(cxx2+1:end)];
    %cc2 = [female_x(1:cxx1) male_x(cxx1+1:cxx2) female_x(cxx2+1:end)];
    
    % Perform single point cross-over
    % cc1,cc2 are bits between xx1 and x22
    cc1 = [male_x(1:cxx_one) female_x(cxx_one+1:end)];
    cc2 = [female_x(1:cxx_one) male_x(cxx_one+1:end)];
    
    % ----------------------- Crossover operators ----------------
    
    % Perform bit-flip mutation on cc1 and cc2
    % Turn this part on for mutation
    %fprintf("Applying mutation operator ------\n")
    cc1 = bit_flip_mutation(cc1);
    cc2 = bit_flip_mutation(cc2);
    
    % Decode scheme
    % Combine all parts to Cn*bit_count array of bits
    child1 = [p1,cc1,p3];
    child2 = [p2,cc2,p4];
    
    % Read b_count number of bits at a time and find decimal value
    % HARDCODED, works only for two children
    child1 = BintoDec(child1, bit_count, look_up);
    child2 = BintoDec(child2, bit_count, look_up);
    
end
% --------------------------- Main genetic function ----------------------

function [p1,p3] = preserve_pos(loc_bunpack, look_up)
    % Force starting position
    p1 = loc_bunpack(1,1);
    p1 = convert_to_bin(p1,look_up);
    % Force ending position
    p3 = loc_bunpack(1,2);
    p3 = convert_to_bin(p3,look_up);
end

% ------------------------------------------------------------------------
function [look_up] = load_lookup_table(bx,bit_count)
look_up = zeros(bx,bit_count); % Initialize, bx by bit_count matrix
point_list = linspace(1,bx,bx); % Row vector
%HARDCODED, but must have a way to generate
    for ii = 1:bx
        var_dec = point_list(1,ii);
        % Before 2019 version
        % bin_arr = de2bi(var_dec,bit_count,'left-msb');
        % Fix for MATLAB 2020 and above.
        % Requires Data Acquisition and Statistics and Machine Learning
        % toolbox
        bin_arr = decimalToBinaryVector(var_dec,bit_count,'MSBFirst');
        
        % Experimental fix by deepseekv3
        %bin_arr = int2bit(var_dec, bit_count); % Convert to bits (column vector)
        %bin_arr = bin_arr'; % Transpose to get a row vector
        %bin_arr = flip(bin_arr, 2); % Flip to get 'left-msb' alignment

        % bin_arr = convert_to_bin(var_dec, look_up);
        
        look_up(ii,:) = bin_arr;
    end
look_up = [transpose(point_list), look_up];    
end

% ------------------------------------------------------------------------
function [out_mmx] = DectoBin(in_mmx, bit_count, look_up)
% Function which converts 1 * n_features into 1 by (bit_count * n_feature * 2)
% in_mmx is a row vector

% How many columns i.e. location index to convert?
row_num = size(in_mmx,1);
col_num = size(in_mmx,2);

% Initalize a matrix but don't declare shape
out_mmx = [];

% Cycle through each integer point and calculate binary representation
    for ii = 1:col_num
        var_dec = in_mmx(1,ii);
        %bin_arr = de2bi(var_dec,bit_count,'left-msb');
        %bin_arr = convert_to_bin(var_dec, look_up);
        bin_arr = decimalToBinaryVector(var_dec,bit_count,'LSBFirst');
        out_mmx = [out_mmx bin_arr];
    end
    
%Sanity check
%out_mmx
%size(out_mmx,2)
    
end

% ---------------------------- Helper functions --------------------------
function [out_mmx] = BintoDec(in_mmx, bit_count, look_up)
% Function to convert 'children' path from bits back to decimal
% in_mmx is a row vector
% How many columns i.e. location index to convert?
row_num = size(in_mmx,1);
col_num = size(in_mmx,2);

% Initalize a matrix but don't declare shape
out_mmx = [];

% How many times we will loop to convert all bits to decimal?
cyc_count = col_num/bit_count; % This must be an integer and multiple of 2

% Global bit position array
% Load base array
glob_count = linspace(1,bit_count,bit_count); % i.e @4bit, arr-> 1,2,3,4  
px = zeros(1,bit_count); % Initialize chunk array

% If the number does not match here, we will set it to 1
look_first_col = look_up(:,1);

% Loop to cycle through each gene, currently in bits of 'bit_count' size
for ii = 1:cyc_count
    if (ii == 1)
        glob_count = glob_count + 0; % This is for first four bits only
    else
        glob_count = glob_count + bit_count; % Update position array
    end
    % Assign glob_count position bits to px
    px(1,:) = in_mmx(1,glob_count);
    % Convert array px back to decimal
    bin_to_var = bi2de(px,'left-msb');
    
    res = ismember(bin_to_var, look_first_col);
    
    if (res == 1)
        % Push this decimal back to output array
        out_mmx = [out_mmx, bin_to_var];
    else
        bin_to_var = 1; % Forcefully send binary bits for 1
        out_mmx = [out_mmx, bin_to_var];
    end
end
%out_mmx % Sanity check
end

% ------------------------------------------------------------------------
function[bit_out] = convert_to_bin(dec_in, look_up)
    % HARDCODED
    look_first_col = look_up(:,1);
    look_bits = look_up(:,[2:end]);
    % Test if dec_in is in this array
    res = ismember(dec_in, look_first_col);
    if (res == 1)
        % Number is found, convert to binary
        dec_loc = find(look_first_col == dec_in);
        bit_out = look_bits(dec_loc,:);
    else
        bit_out = look_bits(1,:); % Forcefully send binary bits for 1
    end
end

% ------------------------------------------------------------------------
function [cc_mutated] = bit_flip_mutation(cc_x)
% Function to perform bit flip mutation
    % Note, cc_x is a row vector
    % Number of bits
    n_bits = size(cc_x,2);

    % Define ranges, HARDCODED, best practice values
    min_pm = 0.001;
    %max_pm = 0.006; % 0.6%
    %max_pm = 0.060; % 6%
    max_pm = 0.600; % 60%
    
    % Generate nn_chromo mutation probabilities between min_pm to max_pm
    prob_mutation = random_generator(n_bits,min_pm,max_pm);
    
    % Shuffle to increase chance of randomization
    prob_mutation = randomize_array(prob_mutation);
    
    % Pick one probability randomly
    rng shuffle; % Reset random seed generator
    prob_mutation = randsample(prob_mutation, 1, 'true', prob_mutation);
    
    % Create a random probability for each bit location
    bit_mutation_prob = random_generator(n_bits,min_pm,max_pm);
    
    % Shuffle to increase chance of randomization
    bit_mutation_prob = randomize_array(bit_mutation_prob);
    
    % Choose a random bit location
    rng shuffle;
    bit_loc = linspace(1,n_bits,n_bits);
    bit_loc = randsample(bit_loc,1);
    bit_prob = bit_mutation_prob(1,bit_loc);
    
    % HARDCODED for two variable problem
    % Is bit_prob greater than prob_mutation? Yes, flip that bit
    if (bit_prob > prob_mutation)
        %fprintf("Mutation performed\n");
        cc_x(1,bit_loc) = ~cc_x(1,bit_loc);
    else
        %fprintf("No mutation!\n");
    end 
    
    % Return mutated array
    cc_mutated = cc_x;
end
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
function [out_vector] = randomize_array(in_vector)
% Function to randomize location of elements in a row vector
% Test if in_vector is a row vector
V = isrow(in_vector);
    if (V == 0)
        in_vector = transpose(in_vector); % Make sure the in vector is a row vector
    end
rand_pos = randperm(length(in_vector)); %array of random positions
% new array with original data randomly distributed 
    for k = 1:length(in_vector)
        out_vector(k) = in_vector(rand_pos(k));
    end
% If our priginal vector was a column vector, we need to put it back in
% correct shape
    if (V == 0)
        out_vector = transpose(out_vector); % Revert it back
    end
end

% ---------------------------- Helper functions --------------------------









