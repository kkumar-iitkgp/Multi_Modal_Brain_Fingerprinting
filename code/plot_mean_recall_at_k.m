function plot_mean_recall_at_k(data_mat,x_array,figure_title,cell_legend,leg_location,cell_label,save_plot_name)


%%
% Summary:
%         1. MATLAB script to generate mean recall @ k plots
%         2. We have kept x-axis-labels for upto 50 (hard coded); needs to
%         be modified for a different range
%
%%
% Function Parameters:
%         Input:
%               1. data_mat: each row contains mean recall@k data 
%               2. x_array: range of values of k 
%               3. cell_legend: legend array
%               4. leg_location: legend location parameter
%               5. cell_label: labels for x-axis and y-axis respectively
%               6. save_plot_name: filename (with path) for saving the
%               plot, if empty: plot won't be saved
%         Output:
%               1. Plot displayed and saved (if filename is not empty)
%
%%
% Reference: 
%           Multi-modal brain fingerprinting: a manifold approximation based framework
% Authors: 
%          Kuldeep Kumar (kkumar@livia.etsmtl.ca), 
%          Laurent Chauvin
%          Matthew Toews (Matthew.Toews@etsmtl.ca) 
%          Olivier Colliot and 
%          Christian Desrosiers (christian.desrosiers@etsmtl.ca)
%     
% LIVIA, ETS Montreal, Canada
% January 2018
%
%%
    Line_width = 3 ;
    Marker_size= 14;
    leg_FontSize=13;
    gca_FontSize = 28;

    color_option{1} = 'r' ;
    color_option{2} = 'g' ;
    color_option{3} = 'b' ;
    color_option{4} = 'c' ;
    color_option{5} = 'm' ;
    color_option{6} = 'y' ;
    color_option{7} = 'k' ;
    
    marker_option{1} = '*' ;
    marker_option{2} = 'o' ;
    marker_option{3} = 'd' ;
    marker_option{4} = 's' ;
    marker_option{5} = 'x' ;
    marker_option{6} = '^' ;
    marker_option{7} = '+' ;

%%
    Num_curves = size(data_mat,1);

    h_fig = figure;
    title(figure_title);
    
    hold on;
    for i=1:Num_curves
        plot(x_array,data_mat(i,:),strcat('-',color_option{i},marker_option{i}),'LineWidth',Line_width,'MarkerSize',Marker_size);
    end

    hold off;
    box on;
    leg_1 = legend(cell_legend,'Location',leg_location);

    set(leg_1,'FontSize',leg_FontSize);
   
    flag_y_lim = 0;  % adaptive mean recall value range for better visualization
    if(flag_y_lim==1)
        ylim([0 1]);
    end
    xlim([0 51]);
    set(gca,'xtick',[0 10 20 30 40 50]);
    xlabel(cell_label{1});
    ylabel(cell_label{2});
    
    set(gca,'FontSize',gca_FontSize);

    set(h_fig,'Position',[50,50,1000,625]);

    if(~isempty(save_plot_name))
        req_rez =1500;
        print(h_fig,save_plot_name,'-dpdf',['-r',num2str(req_rez)],'-opengl');
    end


end