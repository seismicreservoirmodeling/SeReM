function [] = generate_histograms(logs, coloraxis,names)
    figure;
    plot_size_x = size(logs,2);
    ind = 0;
    for i = 1:1:plot_size_x
        for j = 1:1:plot_size_x
            ind = ind + 1;
            subplot(plot_size_x,plot_size_x,ind);
            if i == j
                histogram(logs(:,i) , 'EdgeAlpha',0) ;
            else
                histogram2(logs(:,j),logs(:,i),'FaceColor','flat' , 'EdgeAlpha',0) ;
                view([0 90])
                if nargin>1
                    caxis(coloraxis)
                end
            end            
            if  i == plot_size_x
                if nargin>2
                    xlabel(names{j}, 'FontSize', 12)
                else
                    xlabel({'Z_', j})
                end
            end
            if j == 1
                if nargin>2
                    ylabel(names{i}, 'FontSize', 12)
                else
                    ylabel({'Z_', i})
                end
            end
        end
    end
end
