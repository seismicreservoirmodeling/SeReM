function [x_est, x_var] = Kriging_options(xcoord, dcoords, dvalues, xvar, l, type, krig, angles)
%%  This function works only by limiting the number of max points for conditioning sorted by distande.
% TO DO: % Kriging_options is being implemented to account for additional options
% such as max number of conditioning points, searching neighborhood and type
% of kriging, including the for loop over all the points 

x_est = zeros(size(xcoord,1),1);
x_var = zeros(size(xcoord,1),1);

nmax = 30;
xmean = mean(dvalues);
if size(dvalues,1) > nmax
    
    for i=1:size(xcoord,1)        
        dist = sqrt(sum((xcoord(i,:) - dcoords).^2,2));
        [~, sorted_indices] = sort(dist);
        lowest_indices = sorted_indices(1:nmax);
        if krig == 0
            [x_est(i,:), x_var(i,:) ] = SimpleKriging(xcoord(i,:), dcoords(lowest_indices,:), dvalues(lowest_indices,:), xmean, xvar, l, type, angles);
        else
            [x_est(i,:), x_var(i,:) ] = OrdinaryKriging(xcoord(i,:), dcoords(lowest_indices,:), dvalues(lowest_indices,:), xvar, l, type, angles);
        end
    end
    
else
    
    for i=1:size(xcoord,1)
        if krig == 0
            [x_est(i,:), x_var(i,:) ] = SimpleKriging(xcoord(i,:), dcoords, dvalues, xmean, xvar, l, type, angles);
        else
            [x_est(i,:), x_var(i,:) ] = OrdinaryKriging(xcoord(i,:), dcoords, dvalues, xvar, l, type, angles);
        end
    end
    
end


end