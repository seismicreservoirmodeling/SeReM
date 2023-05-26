function [] = plot_histogram(p, cat)

% Get unique integer values from cat variable
unique_cat = unique(cat);

% Set up colors for plotting
colors = lines(numel(unique_cat));

% Plot histograms for each unique integer value of cat
hold on;
for i = 1:numel(unique_cat)
    current_cat = unique_cat(i);
    current_p = p(cat == current_cat);  % Filter p values for the current cat
    
    histogram(current_p, 'FaceColor', colors(i, :));
end
hold off;

% Add title and labels
title('Histograms of p for different cat values');
ylabel('Frequency');

% Add a legend
legend(cellstr(num2str(unique_cat)), 'Location', 'best');
