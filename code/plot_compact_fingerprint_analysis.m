function plot_compact_fingerprint_analysis(data_mat,num_eig_vec_array,cell_legend,leg_location,cell_label,save_plot_name)


%%
% Summary:
%         1. MATLAB script to generate d-prime/p-val vs num of eig vecs
%            plot for compact fingerprint analysis
%
%%
% Function Parameters:
%         Input:
%               1. data_mat: each row contains mean recall@k data 
%               2. num_eig_vec_array: array for num_eig_vecs for compact
%               fingerprint
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
    Marker_size= 14; %#ok<*NASGU>
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
    hold on;
    for loop_curve = 1:Num_curves         
        plot(data_mat(:,loop_curve),strcat('-',color_option{loop_curve}),'LineWidth',Line_width);
    end
    hold off; 
    box on;
    
       
    set(gca,'xtick',1:length(num_eig_vec_array));
    set(gca,'XTickLabel',num2cell(num_eig_vec_array),'FontSize', 16);
          
    xlabel(cell_label{1,1},'FontSize',gca_FontSize);
    ylabel(cell_label{2,1},'FontSize',gca_FontSize);
     
    leg_1 = legend(cell_legend,'Location',leg_location);
    set(leg_1,'FontSize',leg_FontSize);

 
    set(gca,'FontSize',gca_FontSize);

    set(h_fig,'Position',[50,50,800,500]);

    if(~isempty(save_plot_name))
        req_rez =1500;
        print(h_fig,save_plot_name,'-dpdf',['-r',num2str(req_rez)],'-opengl');
    end


end