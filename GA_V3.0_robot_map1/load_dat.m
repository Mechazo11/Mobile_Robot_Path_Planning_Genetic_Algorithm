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

% loc_xy = [
%     1 7; 
%     1 11; 
%     3 14; 
%     3 1; 
%     5 8; 
%     6 11; 
%     6 4; 
%     8 4; 
%     10 1; 
%     10 7; 
%     10 11; 
%     11 14; 
%     13 12;
%     12 2; 
%     14 8;];

loc_xy = [
    1 7; 
    1 11; 
    3 14; 
    3 1; 
    5 8; 
    6 11; 
    6 4; 
    8 4; 
    10 1; 
    10 7; 
    10 11; 
    11 14; 
    13 12;
    12 2; 
    14 8;];

point_mat = [point_mat_ls, loc_xy];

path_index = [
    1 2 4 5 0;
    2 1 3 0 0;
    3 2 6 12 0;
    4 1 7 9 0;
    5 1 6 7 0;
    6 3 5 11 0;
    7 4 5 8 0;
    8 7 10 14 0;
    9 4 14 0 0;
    10 8 11 15 0;
    11 6 10 12 0;
    12 3 11 13 0;
    13 12 15 0 0;
    14 8 9 15 0;
    15 15 10 13 14;];
end


