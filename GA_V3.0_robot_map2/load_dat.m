function [point_mat, path_index, point_ls, bit_count] = load_dat(start_loc, finish_loc)

% HARDCODED, CHANGES WITH MAP OF THE ENVIRONMENT
bit_count = 4; % [4 = paper, 8 = experimental]
m = 7; % Number of static obstacles, 
start_point = start_loc; % Location in point index
finish_point = finish_loc; % Location in point index

% HARDCODED, CHANGES WITH MAP OF THE ENVIRONMENT
% Range for point index
min_point_index = 1;
max_point_index = 15;

% Return an array of useful points for later use
point_ls = [start_point,finish_point,min_point_index,max_point_index];

% Function to load point_mat and path_index [LATER]
% Generate tables showing point encoding and admissible path
% HARCODED, needs to be adjusted for each static map
point_mat_ls = linspace(1,max_point_index,max_point_index);
point_mat_ls = transpose(point_mat_ls);

loc_xy = [
    1 1;%1 
    3 3;%2 
    6 2;%3 
    9 3;%4 
    6 5;%5 
    9 7;%6 
    12 8;%7 
    6 8;%8 
    9 12;%9 
    3 10;%10 
    6 14;%11 
    14 1;%12 
    3 6;%13
    1 8;%14 
    14 14;];%15
% Coordinates of via point in actual map
point_mat = [point_mat_ls, loc_xy];

% Array which shows admissible points from a location
% HARDCODED
path_index = [
    1 2 3 12 14 0;
    2 3 5 13 14 0;
    3 1 2 4 12 0;
    4 3 5 12 0 0;
    5 2 4 6 8 0;
    6 5 7 8 0 0;
    7 6 9 12 15 0;
    8 5 6 9 10 13;
    9 7 8 11 15 0;
    10 8 11 14 0 0;
    11 9 10 15 0 0;
    12 1 3 4 7 15;
    13 2 8 14 0 0;
    14 1 2 10 13 0;
    15 7 9 11 12 15;];
    % We allow the robot to sit idle in last point


% Note, in finish point we allow robot to linger cause its already at
% finishing point


end


