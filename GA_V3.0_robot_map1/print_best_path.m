function [fig_out] = print_best_path(most_fit_g1, most_fit_value, gen_count, point_mat, point_ls, path_index)
    % HARDCODED, Plot to view best path
    % path_index - table containing path id and admissible points
    % point_mat -- a matrix containing x,y values for each index point
    % point_ls = row vector containing start_point, finish_point, min and
    % https://www.mathworks.com/matlabcentral/answers/156040-how-to-connect-two-points-with-a-line
    
    % Load most fit path
    % Convert most_fit_g1 into a column vector
    cell_path = num2str(most_fit_g1); % For title of plot
    most_fit_g1 = transpose(most_fit_g1); % For other purposes
    
    % Plot the map
    ff1 = figure;
    x_map_cor = point_mat(:,2);
    y_map_cor = point_mat(:,3);
    ff1 = scatter(x_map_cor,y_map_cor,'filled');
    
    fname = append("Generation",' ',num2str(gen_count),' ',"best path");
    val_name = append('Fitness -> ', num2str(most_fit_value));
    xlabel("X direction");
    ylabel("Y direction");
    %title(fname);
    title({fname,cell_path,val_name});

    % Create point labels
    lss = transpose(linspace(1,15,15));
    labels = zeros(15,1);
    for ii = 1:numel(lss)
        labels(ii,1) = string(lss(ii,1));
    end
    hold on;
    
    % Update labels
    ff1 = labelpoints(x_map_cor,y_map_cor,labels,'N',0.15);
    
    % Generate the allowable paths between the points
    
    % How many via points we have
    via_num = size(path_index,1);
    % Load each index
    for pp = 1: via_num
        lss = path_index(pp,:); % Slice out index rpw
        fp = lss(1,1); % Get index value, first position
        % What are the coordinates for fp?
        xx = point_mat(fp, [2,3]);
        
        % Get admissible vertices
        lss_sp = lss(1,[2:end]); 
        % Filter out filler zero
        lss_sp = lss_sp(1,[find(lss_sp ~= 0)]); % Row vector
        
        % Number of lines == number of vertices
        vertex_cnt = size(lss_sp,2);
        
        % Loop through each vertex and update the figure window
        for ii = 1: vertex_cnt
            vxx = lss_sp(1,ii); % Pick a vertex
            yy = point_mat(vxx, [2,3]);
            ff1 = plot([xx(1), yy(1)],[xx(2), yy(2)], '-k');
        end
    end
    
    % How many genes in given path?
    col_cnt = size(most_fit_g1,1); % Remember, this is a column vector now
    % How many segments?
    seg_count = col_cnt - 1;
    
    % Array to hold x and y value of path
    loc_mat = zeros(col_cnt,2); % row x col --> gene_count x 2
    
    % Load x and y values
    for kk = 1: col_cnt
        num_now = most_fit_g1(kk,1);
        loc_mat(kk,:) = point_mat(num_now,[2 3]);
    end
    
    % Plot lines showing path
    xt = loc_mat(:,1);
    yt = loc_mat(:,2);
    ff1 = plot(xt,yt, '-r');
    legend('via points','allowable path','Location','NorthEastOutside','AutoUpdate','off');
    axis auto
    % Plot arrow heads
    for ii = 1:seg_count
        x1 = xt(ii,1);
        y1 = yt(ii,1);
        x2 = xt((ii+1),1);
        y2 = yt((ii+1),1);
        P1=[ x1 y1];
        P2=[ x2 y2];
        D = P2 - P1;
        h = quiver( P1(1), P1(2), D(1), D(2), 0, 'color','red', 'MaxHeadSize',0.5);
        %set(h,'MaxHeadSize',1,'AutoScaleFactor',0.1);
        set(h,'AutoScaleFactor',0.9);
        %legend('via points','allowable path','robot path','Location','NorthEastOutside','AutoUpdate','off');
        % Clear for next set of points
        P1 = [];
        P2 = [];
        D = [];
    end
    
    % Show starting and ending locations
    % Load coordinates for starting and ending position
    s_pos = point_mat(point_ls(1,1),[2 3]);
    f_pos = point_mat(point_ls(1,2),[2 3]);
    % Plot on figure
    ff1 = scatter(s_pos(1,1),s_pos(1,2),100,'og','filled');
    ff1 = scatter(f_pos(1,1),f_pos(1,2),100,'or','filled');
    
    % Increase size of plot
    xlim([0,15]);
    ylim([0,15]);
    
    % Increase size of figure window, change as needed
    x0=10;
    y0=10;
    width=700;
    height=500;
    ff1 = set(gcf,'position',[x0,y0,width,height]);
    
    hold off
    % Push out completed plot
    fig_out = ff1;
end
