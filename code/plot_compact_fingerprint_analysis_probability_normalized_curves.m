function plot_compact_fingerprint_analysis_probability_normalized_curves(cell_data,modality,save_plot_name)


%%
% Summary:
%         1. MATLAB script to compare compact fingerprints for twin/sibling
%         types for a given modality
%         2. Plot probability normalized curves (gamma histogram fitting)
%         for Euclidean distance between twin/sibling pair compact fingerprints.
%
%%
% Function Parameters:
%         Input:
%               1. cell_mat: data cell containing pairwise Euclidean
%               distance between compact fingerprint pairs for MZ, DZ and
%               FS
%               2. modality: name of the modality
%               3. save_plot_name: filename (with path) for saving the
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

    % fitdist and then plot pdf
    h_fig = figure; 
        
    n_bins = 100 ;
    var_dist = 'gamma';
    x_values = 0:0.05:25;

    %FS
    fit_1 = fitdist(cell_data{3,1}(:),var_dist);
    y_1 = pdf(fit_1,x_values);

    %DZ
    fit_2 = fitdist(cell_data{2,1}(:),var_dist);
    y_2 = pdf(fit_2,x_values);

    %MZ
    fit_3 = fitdist(cell_data{1,1}(:),var_dist);
    y_3 = pdf(fit_3,x_values);

    plot(x_values,y_1/1000,strcat('-',color_option{3}),'LineWidth',Line_width);

    hold on
    plot(x_values,y_2/1000,strcat('-',color_option{2}),'LineWidth',Line_width);
    plot(x_values,y_3/1000,strcat('-',color_option{1}),'LineWidth',Line_width);
    hold off;
    
    legend('FS','DZ','MZ','location','northwest');
    set(gca,'FontSize',gca_FontSize);
    ylabel('Probability','FontSize',32);
    xlabel('Euclidean Dist','FontSize',32);
    %ylim([0 0.05]);
    %xlim([4 17]);
    title(['Probability Normalized plot: ' modality{1}]);
    
    set(h_fig,'Position',[50,50,1000,625]);

    if(~isempty(save_plot_name))
        req_rez =1500; print(h_fig,save_plot_name,'-dpdf',['-r',num2str(req_rez)],'-opengl');
    end
    
end       