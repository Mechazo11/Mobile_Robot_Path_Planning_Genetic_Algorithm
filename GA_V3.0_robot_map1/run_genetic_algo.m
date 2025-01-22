function [new_gen_2] = run_genetic_algo(gen_count, bit_count,X1, N, path_index, point_mat, point_ls,m)
    
    % UNIQUE PROBLEM FUNCTION / ALGO FUNCTION
    % ------------------------------------------------------------------------
    % Input
    % gene_count - Integer number for current generation
    % bit_count - number of bits to encode integer values to binary
    % X1 -- Current generation matrix
    % N -- Number of candiates per generation
    % path_index - table containing path id and admissible points
    % point_mat -- a matrix containing x,y values for each index point
    % point_ls = row vector containing start_point, finish_point, min and
    % mix index integers -- IMPORTANT for creating children
    
    % Output
    % new_gen_2 -- matric containing (N+1) candidates
    % -----------------------------------------------------------------------
    
    % We need to know the finish point
    %finish_point = point_ls(1,2); % point_ls is a row vector
    
    % Step 0.1: Survival of the fittest, drop row with least fit from N+1
    % candidates
    % Brings X1 back to N number of rows
    % HARDCODED
    row_in = size(X1, 1);
    if (row_in > N)
        [X1] = bring_back_N(X1, path_index, point_mat);
    end
    
    % Experimental operator: Check if we have more than 80% similarity
    X1 = check_minima_loc(X1,bit_count,N, path_index, point_mat, point_ls,m);
    
    % Step 0.5: Global elite operator
    % Inject the overall best path in the current iteration
    X1 = inject_global_best(X1, path_index, point_mat);
    
    % Need a table to show these answers later
    % DEBUG
    % fprintf("\n---------------------- Candiates going in----------------\n")
    % X1
    % fprintf("\n---------------------- Candiates going in----------------\n")  
    
    % Step 1b: Test if al fit_g1 is zero, UNIQUE to this problem
    % fit_g1 is a column vector, totat_fit_g1 is a scalar
    [fit_g1, total_fit_g1] = eval_obj(X1, path_index, point_mat);
    
    % DEBUG
    fprintf("\n---------------------- Candiates fitness----------------\n")
    fit_g1'
    fprintf("\n---------------------- Candiates fitness----------------\n")
    
    % Step 2: Remember best fit individual from this generation
    % What happens if all values are zero?
    % How many columns before adding fitness values
    [most_fit_g1, fit_xx] = find_most_fit(X1, path_index, point_mat);
    
    % Push this generation's most fit candidate back to global space for
    % plotting
    % I realize using global variables may not be a good idea but this
    % idea works
    global most_fit_candidate;
    global current_candidate_fitness;
    most_fit_candidate = most_fit_g1;
    current_candidate_fitness = fit_xx;
    
    % Step 3 Evaluate fraction of fitness and cumulative probability
    [frac_fit_g1, cumu_prob_g1] = eval_fraction_fitness(fit_g1, total_fit_g1);
    
    % Step 3a and 3b Tabulate results, will be done later
    
    % temporary
    %fprintf("\n---------------------------------------################ -----------------------------------\n");
    %x_test = [X1,fit_g1,frac_fit_g1, cumu_prob_g1]
    %format short;
    %x_test = [X1 fit_g1]; % Temporary 
    %avg_fit_g1 = total_fit_g1 / N;
    %x_test_2 = table(avg_fit_g1)
    %fprintf("\n---------------------------------------################ -----------------------------------\n");
    
    % Step 4 Select N random numbers for numerical wheel of fortune
    rand_ls = random_generator((N*2), 0.25, 0.99);
    comp_rand = datasample(rand_ls, N,'Replace',false);
    
    % Step 5 Determine mate using wheel of fortune
    mate_matrix = find_mates_with_replacement(X1, cumu_prob_g1, comp_rand);
    %mate_matrix = find_mates_without_replacement(X1, cumu_prob_g1, comp_rand);
    
    % Step 6 Reshape to have candiates form couple per row
    mmx = reshape_long_row(mate_matrix);
    
    % Step 7 Perform cross over, massive overhaul
    new_gen_2 = create_new_gen(mmx,bit_count,point_ls);

    % Step 8 Reshape new generation matrix back to N by num_of_gene
    new_gen_2 = reshape_back(new_gen_2);
    
    % Step 9 Perform elitism, copy-paste the strongest guy from the last generation
    new_gen_2 = [new_gen_2;most_fit_g1];
    
end

function [X1_renew] = check_minima_loc(X1,bit_count,N, path_index, point_mat, point_ls,m)
    % Function which reinitializes the candiadate matrix if more than 80%
    % candiate paths are same
    % Default value to pass
    X1_renew = X1;

    % Calculate fitness
    [fit_g1, ~] = eval_obj(X1, path_index, point_mat);
    
    %fit_g1 % Debug purpose
    % Make sure fit_g1 is a column vector
    if (isrow(fit_g1))
        % Transpose to get a column vector
        fit_g1 = transpose(fit_g1);
    end
    
    % 80% of value of N
    chk_xx = int16(N - (0.8 * N));
    
    % How many unique fitness?
    u_fit = numel(unique(fit_g1));
    
    % Check condition
    if (u_fit <= chk_xx)
        % We are locked in a local minima, reshuffel
        [X1,~] = random_g1(N,bit_count,m,point_ls(1,1), point_ls(1,2), point_ls(1,3), point_ls(1,4));
        % Push new update
        X1_renew = X1;
    end
end

% --------------------------------------------------------------------------
function [X1_best_global] = inject_global_best(X1, path_index, point_mat)
    % UNQIUE PROBLEM FUNCTION, HARDCODED
    % Load global variable
    global global_elite_count;
    global global_elite_path;
    global global_elite_fitness
    
    % Load working variables
    gg_current_fit = 0;
    gg_current_path = [];
    
    % Calculate fitness
    [fit_g1, ~] = eval_obj(X1, path_index, point_mat);
    
    %fit_g1 % Debug purpose
    % Make sure fit_g1 is a column vector
    if (isrow(fit_g1))
        % Transpose to get a column vector
        fit_g1 = transpose(fit_g1);
    end
    
    % Update global elite path
    if (global_elite_count == 1)
        % In generation 1, find the best path
        % Which path is most fit?
        [global_elite_fitness,I] = max(fit_g1);
        global_elite_path(1,:) = X1(I,:);
        % Update count number
        global_elite_count = global_elite_count + 1;
    else
        % In generation n,
        % Find current path's elite member
        % Compare if current member's fitness is bigger than global elite's fitness
        [gg_current_fit,ii] = max(fit_g1);
        if(gg_current_fit > global_elite_fitness)
            % Update global fitness
            global_elite_fitness = gg_current_fit;
            % Update global fitness path
            global_elite_path(1,:) = X1(ii,:);
            global_elite_count = global_elite_count + 1;
        end
    end
    
    % Push global best into current candiate array
    %X1
    %fit_g1
    X1(1,:) = global_elite_path(1,:);
    X1_best_global = X1;
end

% --------------------------------------------------------------------------
function [X1_best] = bring_back_N(X1, path_index, point_mat)
    % UNQIUE PROBLEM FUNCTION, HARDCODED
    % Remove weakest design from the candidate space
    
    % Calculate fitness
    [fit_g1, ~] = eval_obj(X1, path_index, point_mat);
    
    %fit_g1 % Debug purpose
    % Make sure fit_g1 is a column vector
    if (isrow(fit_g1))
        % Transpose to get a column vector
        fit_g1 = transpose(fit_g1);
    end
    
    % Which index is minimum value
    [M,I] = min(fit_g1);
    X3_1 = X1(1:(I-1),:);
    X3_2 = X1((I+1):end,:);
    X3 = [X3_1;X3_2];
    
    % Return candidates with N rows
    X1_best = X3; 
end

% ---------------------------------------------------------------------------
function [mmx_reshape] = reshape_back(new_gen)
    % Here new_gen is a (N/2) * (2*num_gene) shape

    % Row and column numbee in new_gen matrix
    rx = size(new_gen,1);
    cx = size(new_gen,2);
    
    % How many rows we need to make, number of col in each row?
    row_x = rx * 2; % Each new_gen row has a set of child
    col_x = cx / 2; % There are 2 times the genes
    
    % Initialize a temporary matrix
    %temp_mat = zeros(row_x, col_x);
    temp_mat = [];
    
    % Array which will control row positions of tempo_mat
    tempo_rx = linspace(1,rx,rx); % Seed array
    iix = 1;
    
    % Cycle through each row of new_gen matrix 
    for ii = 1:rx
        ttx = new_gen(ii,:);
        ttx = reshape_back_one_set(ttx);
        temp_mat = [temp_mat;ttx];
        %temp_mat(tempo_rx,:) = ttx
        %tempo_rx = tempo_rx + 2; % i.e if 1,2 was finished go to 3,4
    end
    % Return reshape matrix
    mmx_reshape = temp_mat;
end

% ---------------------------------------------------------------------------
function [mmx_reshape] = reshape_back_one_set(new_gen)
% Function which reshapes back to N by num_gene shape
% new_gen must be 1 * (2 * num_gene) shape
mmx = new_gen; % Careful, this has to be a row vector
szz1 = size(new_gen,1);
szz2 = size(new_gen,2);
% Define a temporary matrix of shape 2 by num_gene shape
ppx = zeros((szz1*2),(szz2 / 2));

% Record number of rows and columns of the reshape matrix
ppx_row = size(ppx,1);
ppx_col = size(ppx,2);
ctn = 1;
  % This nested loop is actually extracting the genes's for each candidate
  % from the flattened matrix
  for pp_row = 1: ppx_row % Row selector
      for pp_col = 1:ppx_col % column selector
          ppx(pp_row, pp_col) = mmx(1,ctn); 
          ctn = ctn + 1;
       end
  end
mmx_reshape = ppx;
end

% ---------------------------------------------------------------------------
function [fitness, total_fit] = eval_obj(X, path_index, point_mat) 
    % UNIQUE PROBLEM FUNCTION
    
    % Determine shape of matrix X
    row_X = size(X,1);
    col_X = size(X,2);
    
    % 2D localization
    coordinate_num = 2;
    
    % Initialize a colum vectors to hold fitness values and cumu_fit
    fitness = zeros(row_X,1);
    total_fit = zeros(row_X,1);
    
    % Calculate fitness for each design path
    for ii = 1:row_X
        % Load up the design
        %fprintf("\nCurrent design idx --> %d", ii);
        %X(ii,:) % Debug
        current_path = X(ii,:);
        fitness(ii,1) = eval_path(current_path, path_index, point_mat);   
    end
    
    % Multiply by 1000 to make it integer
    fitness = fitness .* 1000; % Fitness here is column vector
    
    % Add up all fitness to find total fitness
    total_fit = sum(fitness,1); % fitness is a column vector
end

% ---------------------------------------------------------------------
function [fitness] = eval_path(current_path, path_index, point_mat)
% UNIQUE PROBLEM FUNCTION
row_X = size(current_path,1);
col_X = size(current_path,2);

% How many line segments?
seg_count = col_X - 1;

% Create a row vector to store individual sums
seg_sum = zeros(1,seg_count);

% Extract first column from path_index and convert it to array
point_idx = path_index(:,1);

% HARDCODED
% Create a 2 by 2 square matrix to hold X and Y values for two points
loc_mat = zeros(2,2); 

for ii = 1: seg_count
     % Load two points in sequence
     %fprintf("\nLine segment --> %d\n",ii);
     fp = current_path(1, ii);
     sp = current_path(1, (ii+1));
     % View values
     % Load up the admissible point array for fp
     admit_arr = path_index(fp, 2:end);
     
     % is sp a valid vertex?
     is_valid  = valid_vertex(admit_arr,sp);
     
     if (is_valid == 1)
        % Fill in loc_mat with coordinate values
        % HARDCODED 2nd and 3rd columns are x and y values respectively
        loc_mat(1,:) = point_mat(fp, [2 3]);
        loc_mat(2,:) = point_mat(sp, [2 3]);
        %fprintf("\nfp --> %d, sp--> %d\n", fp,sp)
        %loc_mat % Debug
        seg_sum(1,ii) = sqrt((loc_mat(2,1) - loc_mat(1,1))^2 + (loc_mat(1,2) - loc_mat(2,2))^2);
        
        % Penal factor
        pp_fact = 0;% Experimental number, 0 no penalty
        
        % Back tracking, we penalize the path
        if (fp > sp)
            % Multiply current sum with (-1)
            seg_sum(1,ii) = seg_sum(1,ii) + 0;    
        end
        
        % go to next line segment
        continue;
     else
         % We don't have any more valid line segments
         % Clean out all sums cuz the path is infesible
         seg_sum = 1000; % Experimental
         break;
     end
end

% What is the seg_sum look like now?
p = sum(seg_sum,2);

% Here we have all the distance values for each sum
d_tot = sum(seg_sum,2); % sum(A,2) calculates sum along row
fitness = 1/d_tot;

% is fitness 0, Yes, assign value of 0.001, avoid NaN problem
    if(fitness == 0)
       fitness = 0.001;
    end
end 

% ---------------------------------------------------------------------
function [response] = valid_vertex(valid_arr, sp)
    % valid_arr contains the list of allowed locations
    % sp is the point we want to test
    ls_cnt = numel(valid_arr); % How many elements
    chk = 0; % Assume no match
    % Check sp against all allowable vertex
    for ii = 1: ls_cnt
        test_val = valid_arr(1,ii);
        if (sp == test_val)
            chk=1;
            break; % As soon as we get a match, exit loop
        end
    end
    % Update response
    response = chk;
 end
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
function [X_tested] = all_check_stationary(X,path_index,finish_point)
    % UNIQUE PROBLEM FUNCTION
    % Backup orginial matrix
    X_back = X;
    
    % Count number of rows
    row_X = size(X,1);
    
    for ii = 1:row_X
    % Load one path at a time
    current_path = X(ii,:);
    
    % Test if the robot is sitting in same point
    tested_path = check_stationary(current_path, path_index, finish_point);
    
    % Update tested_path into main design space
    X(ii,:) = tested_path;
    end
    
    % Push updated design into X_tested
    X_tested = X;
end

% ---------------------------------------------------------------------
function [tested_path] = check_stationary(current_path, path_index, target_point)
% HARDCODED,
% Idea is if we see an infesible point from a chosen point
% there is 50-50 chance that code will change the point
% to one of the allowable points
cc_old = current_path;
cg_prob = 0.75; % Experimental

% Extract first column from path_index and convert it to array
% Number of index points 
point_idx = path_index(:,1);

% How many genes in one chromosome?
col_X = size(current_path,2);

% How many line segments?
seg_count = col_X - 1;

% Create a large random probability pool
rand_ls = random_generator((4*col_X),0,0.98);

% Load points
    for ii= 1:seg_count
    %fprintf("Seg num --> %d\n", ii);
    fp = current_path(1, ii);
    sp = current_path(1, (ii+1));
    
    % Is fp already at ending point
    if (fp == target_point)
        break; % No need to perform random shuffle
    end
    % Comment the forced chance_prob value, only for debug
    %cc_prob = 0.1;
    cc_prob = randsample(rand_ls, 1, 'true', rand_ls);
    
    % Are they identical
    if(fp == sp)
        %fprintf("\nRobot sitting idle, random swapping\n")
        % Is chance_prob < cg_prob ? Yes change, No go to next line seg
        if (cc_prob < cg_prob)
            % Look up admissible from path index
            sp_who = find(point_idx == sp);
            admit_arr = path_index(fp, 2:end); % Row vector
            % Which admissible vertices are greater than fp?
            admit_arr = find_bigger_vertex(admit_arr, fp);
            %sp_new = randsample(admit_arr, 1, 'true', admit_arr);
            
            sp_new = max(admit_arr);
            % If NaN is chosen, set sp = sp
            res = isnan(sp_new);
            if (res == 1)
                sp_new = sp; % We got a nan, so don't change anything
            else
                sp = sp_new; % Update sp with new value
            end
        end
    end
    % Update points in the current path
    %sp % Why error
        try
            current_path(1, (ii+1)) = sp;
        catch
            fprintf("Catch statement executed");
        end
    end

% Push update
tested_path = current_path;
end

% ---------------------------------------------------------------------
function [big_vertex_arr] = find_bigger_vertex(arr_in, fp)
% UNIQUE PROBLEM FUNCTION
% arr_in --> Array of admissible points for point 'fp'
% fp --> First point for line segment in consideration
% big_vertex_arr --> Array containing vertexes in the forward direction
logical_idx = arr_in(:)>fp;
big_vertex_arr = arr_in(logical_idx);
end

% ---------------------------------------------------------------------
function [feasible] = check_feasible_path(current_path, path_index, loc_bunpack)
 % UNIQUE PROBLEM FUNCTION
 % Define variables to stop MATLAB error
 feasible = zeros(1,1); % This is a scalar term
 
 % Extract first column from path_index and convert it to array
 % All rows, first column 
 point_idx = path_index(:,1); 
 
 % How many line segments?
 col_X = size(current_path,2);
 seg_count = col_X - 1;
 
 % elements of the path_index's 2nd column
 admit_arr = []; % Initiate
 ax = loc_bunpack(1,1);
 bx = loc_bunpack(1,2);
 spx = 0; % Special flag against start-end position issue
 
 for ii = 1: seg_count
     % Load two points in sequence
     %fprintf("\nLine segment --> %d\n",ii);
     fp = current_path(1, ii);
     sp = current_path(1, (ii+1));
     
     % What is the location of first point in path_idx?
     [r_x] = find(point_idx == fp);
     
     % Load up the admissible point array for r_x
     %fprintf("\nLoading admissible array for %d\n",r_x);
     % Uncomment for debug only
     %fp
     %sp
     admit_arr = path_index(fp, 2:end); % Row vector
     
     % RESUME FROM HERE, make new ismember function
     res = ismember(sp, admit_arr);
     
     % Is sp admissible for current fp?
     if (res ==1)
        %fprintf("Point %d in admissible array",sp);
        % Go to next line segment
        continue 
     else
        %fprintf("Point %d not admissible array",sp); 
        feasible = 0;
        spx = 1;
        break % We don't check this path anymore
    end        
 end
 
    if (ii == seg_count)
        % We checked all line segment and found feasible
        %fprintf("\nAll line segment feasible\n");
        feasible = 1;
    end
    
% What is the verdict? Comment for debug only
%feasible
end

% ---------------------------------------------------------------------------
function [X1_best, X1_fitness] = find_most_fit(X1, path_index, point_mat) 
    % Unqiue problem function
    
    % Define variable
    best_mat = [];
    
    % Calculate fitness
    [fit_g1, ~] = eval_obj(X1, path_index, point_mat);
    
    % Which index is maximum value
    [M,I] = max(fit_g1);
    
    % Return value of most fitness
    X1_fitness = M;
    
    % Extract that row
    best_mat = X1(I,:);
    
    % Return path with best fitness
    X1_best = best_mat;
end

% ---------------------------------------------------------------------------
function [frac_to_fit, cumu_prob] = eval_fraction_fitness(fitness, cumu_fit)
    frac_to_fit = fitness./cumu_fit;
    cumu_prob = cumsum(frac_to_fit,1);
end

% ---------------------------------------------------------------------------
function X22 = find_mates_with_replacement(X_current, cumu_prob, comp_rand)
    % Allow a candiate to be repeated in the mating pool
    
    % Find number of columns in design matrix
    X1_col_last_num = size(X_current, 2);
    
    % Make sure cumu_prob and comp_rand are row vectors
    if( isrow(cumu_prob)== 0)
    cumu_prob = transpose(cumu_prob);
    end
    
    if(isrow(comp_rand)==0)
    comp_rand = transpose(comp_rand);
    end
    
    % Delete this later
    %comp_rand = [0.184, 0.1383, 0.7099, 0.7043, 0.0627, 0.7647];
    
    % Put together X1, cumu_prob and comp_rand column wise
    X11 = [X_current, (cumu_prob).', (comp_rand).']; % Big matrix
    NX = size(X11,1); % Number of rows i.e number of candiates
    X22 = zeros(NX,size(X_current,2)); % X22 size must be N by n_features
    num_to_iter = NX;
    
    % HARDCODED, location of cumulative probability and 
    % random probability locations
    cumu_prob_col = X1_col_last_num + 1;
    rand_prob_col = X1_col_last_num + 2;
    
    % Scratch pad variable
    row_prob = 0;
    rand_prob = 0;
    test_idx = 1;
    cnt = 1;
    
    %X11 % Debug
    
    % Main loop
    for ii = 1:num_to_iter
        rand_prob = X11(ii,rand_prob_col); % Choose random number sequentially
        %fprintf("Rand prob now --> %f\n",rand_prob);
        
        % Look through each cumulative probability
        for kk = 1: num_to_iter
            % Load each cumulative probability sequentially
            row_prob = X11(kk,cumu_prob_col);
            %fprintf("Cumulative prob now --> %f\n",row_prob);
            
            % Check if current cumulative prob is bigger than current row
            % vector
            if (row_prob > rand_prob)
                %X22(ii,:) = X_current(kk,:);
                X22(ii,[1:X1_col_last_num]) = X_current(kk,[1:X1_col_last_num]);
                
                % Reset for next random number
                rand_prob = 0;
                row_prob = 0;
                cnt = 1;
                break; % We got a hit
            else
                cnt = cnt + 1; % Dummy to prevent error
                if (cnt<num_to_iter)
                    continue % Go to next random probability value
                else
                    % Guard against overflow
                    % change random probability value
                    % Exceeded num_iter count
                    % Reset cnt
                    % reset kk
                    % continue again
                    rand_prob = 0.111;
                    cnt = 1;
                    kk = 1;
                    continue
                end
            end
        end
    end
    %X22 % What is X22 now?
end

% ---------------------------------------------------------------------------
function X22 = find_mates_without_replacement(X_current, cumu_prob, comp_rand)
    % Find number of columns in design matrix
    X1_col_last_num = size(X_current, 2);
    
    % Make sure cumu_prob and comp_rand are row vectors
    if( isrow(cumu_prob)== 0)
    cumu_prob = transpose(cumu_prob);
    end
    
    if(isrow(comp_rand)==0)
    comp_rand = transpose(comp_rand);
    end
    
    % Put together X1, cumu_prob and comp_rand column wise
    X11 = [X_current, (cumu_prob).', (comp_rand).']; % Big matrix
    NX = size(X11,1); % Number of rows i.e number of candiates
    X22 = zeros(NX,size(X_current,2)); % X22 size must be N by n_features
    num_to_iter = NX;
    
    % HARDCODED, location of cumulative probability and 
    % random probability locations
    cumu_prob_col = X1_col_last_num + 1;
    rand_prob_col = X1_col_last_num + 2;
    
    % Scratch pad variable
    row_prob = 0;
    rand_prob = 0;
    test_idx = 1;
    cnt = 1;
    
    % Counter to keep track of how many counts we did
    for ii = 1: num_to_iter
        row_prob = X11(ii,cumu_prob_col); % Choose the cumu_probability for this candidate
        %fprintf("\nCumu prob now --> %f\n",row_prob);
        %fprintf("------------------------\n");
        
        % Look thru each "random" probability sequentially
        for kk = 1: num_to_iter
            %fprintf("\nCumu prob now --> %f\n",row_prob);
            rand_prob = X11(kk,rand_prob_col); % Get probability for position at kk
            %fprintf("Random prob now --> %f\n",rand_prob);
            %%%Test if current cumu probability is greater than current value of random probability
            if (row_prob >= rand_prob)
                X22(kk,[1:X1_col_last_num]) = X_current(ii,[1:X1_col_last_num]);
                X11(kk,rand_prob_col) = 1.1; % This will prevent this position random number to evaluated again in next round
                %X11(:,rand_prob_col) % Did it actually set to 1.1?
                
                %fprintf("-----------------------------------------------------------\n");
                %fprintf("Assigned G_%d position %d going to next candidate\n", ii,kk);
                %fprintf("-----------------------------------------------------------\n");
                break % We stop the inner loop and move onto next candidate
            else
                cnt = cnt + 1; % Dummy to prevent error
                if (cnt<num_to_iter)
                    continue % Go to next random probability value
                else
                    % change actual probability value of this candidate to 1
                    % Exceeded num_iter count
                    % Reset cnt
                    % reset kk
                    % continue again
                    %fprintf("Boosted current row probability\n")
                    row_prob = 0.99999;
                    cnt = 1;
                    kk = 1;
                    continue
                end
            end
        end
    end
% X22 now is in N x n_feature form. We will reshape them into (N/2,
% n_feature *2) for easier computation of next generation design candidates
end

% ---------------------------------------------------------------------------
function mxx_reshape = reshape_long_row(mate_matrix)
% Function which reshapes N by n_feature to (N/2) by (n_feature*2) matrix
assert(mod(size(mate_matrix,1), 2) == 0, 'Input matrix must have even number of rows');
mmx = mate_matrix;
mmx = mmx'; % Shape --> (n_feature by N)
mmx = reshape(mmx,1,[]); % All elements flattened into a row vector
num_count = size(mmx,2);
szz1 = size(mate_matrix,1);
szz2 = size(mate_matrix,2);
ppx = zeros((szz1/2),(szz2 * 2)); % Temporary array to hold (N/2) by (n_feature*2)

% Record number of rows and columns of the reshape matrix
ppx_row = size(ppx,1);
ppx_col = size(ppx,2);
ctn = 1;
  % This nested loop is actually extracting the genes's for each candidate
  % from the flattened matrix
  for pp_row = 1: ppx_row % Row selector
      for pp_col = 1:ppx_col % column selector
          ppx(pp_row, pp_col) = mmx(1,ctn); 
          ctn = ctn + 1;
       end
  end
mxx_reshape = ppx;
end

