function [ variable_nonParametric ] = uniform_to_nonParametric( data2transform, reference_variables, grid_size, cond_variables)

    n_cond_vars = size(cond_variables,2);
    n_vars = size(reference_variables,2);
    
    min2norm = min(reference_variables);
    
    reference_variables = reference_variables - repmat(min2norm, size(reference_variables,1), 1);
    max2norm = max(reference_variables);
    reference_variables = reference_variables ./ repmat(max2norm, size(reference_variables,1), 1);
    cond_variables = cond_variables - repmat(min2norm(1:n_cond_vars), size(cond_variables,1), 1);
    cond_variables = cond_variables ./ repmat(max2norm(1:n_cond_vars), size(cond_variables,1), 1);
    
    %data2transform = data2transform - repmat(min2norm, size(data2transform,1), 1);
    %data2transform = data2transform ./ repmat(max2norm, size(data2transform,1), 1);
    
    % tic
    %num_point_without_statistic = 0;
    variable_nonParametric = zeros(size(data2transform));
    parfor i = 1:size(data2transform,1) % for each point of the grid, repeat
        reference_variables_filtered = reference_variables;
        data = data2transform(i,:);
        if isnan(cond_variables(i,1))
            
        else
            for j =1:n_vars % for each conditional variable, repeat
                if j <= n_cond_vars
                    variable_nonParametric(i,j) = cond_variables(i,j);
                    index = abs(reference_variables_filtered(:,j) - cond_variables(i,j)) < grid_size; %filtra apenas valores em torno do valor amostrado
                    reference_variables_filtered(~index,:) = NaN;
                else
                
                    % inverse cumulative from the conditional
                    invcumhist = sort(reference_variables_filtered(~isnan(reference_variables_filtered(:,1)),j));  %gera cumulativa
                    
                    if size(invcumhist,1) > 1
                        invcumhist = invcumhist + [1:length(invcumhist)]'*1e-10;  % Sum an infinitesinal line to avoid intep problems
                        
                        domain = [1:length(invcumhist)]/length(invcumhist);
                        variable_nonParametric(i,j) = interp1( domain, invcumhist, data(j) );
                        
                        % Few points to compute the inverse cumulative and the draw value is out side the existent range
                        if isnan(variable_nonParametric(i,j))
                            if data(j)<=min(domain)
                                variable_nonParametric(i,j) = min(invcumhist);
                            else
                                if data(j)>=max(domain)
                                    variable_nonParametric(i,j) = max(invcumhist);
                                end
                            end
                        end
                    else
                        %num_point_without_statistic = num_point_without_statistic + 1;                
                        %disp('Not enough data for contitioning for ' +string(num_point_without_statistic)+' points. The method will draw from the marginal. It might generate artifacts. Consider using KDE to increase the number of data points or increasing the grid_size parameter.')
                        disp('Not enough data for contitioning for specific point. The method will draw from the marginal. It might generate artifacts. Consider using KDE to increase the number of data points or increasing the grid_size parameter.')
                        invcumhist = sort(reference_variables(:,j));  %inverse cumulative from the marginal
                        variable_nonParametric(i,j) = interp1( [1:length(invcumhist)]/length(invcumhist),invcumhist,data(j) );
                    end
                    
                    index = abs(reference_variables_filtered(:,j) - variable_nonParametric(i,j)) < grid_size; %filtra apenas valores em torno do valor amostrado
                    reference_variables_filtered(~index,:) = NaN;
                end
            end
        end
    end
    
    
    variable_nonParametric = variable_nonParametric .* repmat(max2norm, size(variable_nonParametric,1), 1);
    variable_nonParametric = variable_nonParametric + repmat(min2norm, size(variable_nonParametric,1), 1);
    
    
    end
    % toc