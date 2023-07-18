function [ variable_nonParametric ] = uniform_to_nonParametric( data2transform, reference_variables, grid_size)


min2norm = min(reference_variables);

reference_variables = reference_variables - repmat(min2norm, size(reference_variables,1), 1);
max2norm = max(reference_variables);
reference_variables = reference_variables ./ repmat(max2norm, size(reference_variables,1), 1);

%data2transform = data2transform - repmat(min2norm, size(data2transform,1), 1);
%data2transform = data2transform ./ repmat(max2norm, size(data2transform,1), 1);

% tic
num_point_without_statistic = 0;
for i = 1:1:size(data2transform,1) %para cada ponto repete;
    reference_variables_filtered = reference_variables;
    data = data2transform(i,:);
    for j =1:1:size(data,2) %para cada variavel repete
        
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
            num_point_without_statistic = num_point_without_statistic + 1;
            disp('Not enough data for contitioning for ' +string(num_point_without_statistic)+' points. The method will draw from the marginal. It might generate artifacts. Consider using KDE to increase the number of data points or increasing the grid_size parameter.')
            invcumhist = sort(reference_variables(:,j));  %inverse cumulative from the marginal
            variable_nonParametric(i,j) = interp1( [1:length(invcumhist)]/length(invcumhist),invcumhist,data(j) );
        end
       
        index = abs(reference_variables_filtered(:,j) - variable_nonParametric(i,j)) < grid_size; %filtra apenas valores em torno do valor amostrado                        
        reference_variables_filtered(~index,:) = NaN;
    end
end
% toc

variable_nonParametric = variable_nonParametric .* repmat(max2norm, size(variable_nonParametric,1), 1);
variable_nonParametric = variable_nonParametric + repmat(min2norm, size(variable_nonParametric,1), 1);


end