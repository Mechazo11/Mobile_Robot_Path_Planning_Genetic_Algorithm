function [gene_ls, chromo_len] = random_g1(N,bit_count,m, s_pos, f_pos, min_index, max_index)
% UNIQUE PROBLEM FUNCTION
% ------------------------------------------------------------------------
% Input
% N-- number of candidates to genereate
% bit_count -- number of bits per gene
% m -- number of static obstacles
% s_pos -- Starting position index
% f_pos -- Ending position index
% min_index -- integer "1"
% max_index -- furthest point from position "1"

% Output
% rand_ls -- N number of candiates containing feasbile and infeasbile path

% Important variables
% Cn -- number of genes in a candidate

% ------------------------------------------------------------------------

% How many paths per candidate? Find chromosome lenght
[Cn, chromo_len] = chromosome_length(bit_count,m);

% Generate N candidates as a (N * Cn), values in point index(integer)
gene_ls = seed_path(N,Cn,min_index, max_index); 

% Change 1st value correspoinding to chosen start point
    for ii= 1:N
        gene_ls(ii,1) = s_pos;
    end
% Change last value correspoinding to chosen finish point
    for ii= 1:N
        gene_ls(ii,end) = f_pos;
    end
%gene_ls %Debug    

% HARDCODED, inserting two feasible path
% record current number of rows
%1 4 9 14 8 10 15 13 15; % 38.51
fes_1 = [
    1 5 7 4 9 14 8 10 15; % 
    1 5 7 4 9 14 15 10 15; % 27.52
    1 4 9 14 15 13 12 13 15; % 27.77
    ];

% Push the three feasible solutions into the initial mate matrix
% HARDCODED, CHANGES WITH MAP OF THE ENVIRONMENT
gene_ls((end-2),:) = fes_1(3,:);
gene_ls((end-1),:) = fes_1(2,:);
gene_ls((end),:) = fes_1(1,:);
end

% Helper functions
%-------------------------------------------------------
function [Cn, chromo_len] = chromosome_length(bit_count, m)
% UNIQUE PROBLEM FUNCTION
% ------------------------------------------------------------------------
% Input
% bit_count -- number of bits per gene
% m -- number of static obstacles

% Output 
% Cn -- number of genes in a candidate
% chromo_len -- no. of bits in candidates
% ------------------------------------------------------------------------
Cn = m + 2;
chromo_len = int16(Cn * bit_count);
end

%-------------------------------------------------------
function [gene_ls] = seed_path(N,Cn,min_index, max_index)
% UNIQUE PROBLEM FUNCTION
% ------------------------------------------------------------------------
% Input

% Output

% Important variables
% ------------------------------------------------------------------------

gene_ls = zeros(N,Cn); % Initialize
% Row selector
    for i= 1:N
      % Column selector
      for k = 1:Cn
        gene_ls(i,k) = int16(random_generator(1,min_index, max_index));
      end
    end
end

