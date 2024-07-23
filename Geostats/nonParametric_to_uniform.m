function [ variable_uniform ] = nonParametric_to_uniform( data2transform, reference_variables, grid_size)

min2norm = min(reference_variables);

reference_variables = reference_variables - repmat(min2norm, size(reference_variables,1), 1);
max2norm = max(reference_variables);
reference_variables = reference_variables ./ repmat(max2norm, size(reference_variables,1), 1);

data2transform = data2transform - repmat(min2norm, size(data2transform,1), 1);
data2transform = data2transform ./ repmat(max2norm, size(data2transform,1), 1);

n_cond_vars = size(data2transform,2);
% tic
parfor i = 1:1:size(data2transform,1) %para cada simulacao repete;
    reference_variables_filtered = reference_variables;
    nois = data2transform(i,:);
    for j =1:1:n_cond_vars %para cada variavel repete
        
        logs_sub = sort(reference_variables_filtered(~isnan(reference_variables_filtered(:,1)),j));  %gera cumulativa
        
        if size(logs_sub,1) > 1
            logs_cumhist = logs_sub;
            logs_sub2 = nois(j);
            % Sum an infinitesinal line to avoid intep problems
            logs_cumhist = logs_cumhist + [1:length(logs_cumhist)]'*1e-10;
            variable_uniform(i,j) = interp1( logs_cumhist,[1:length(logs_cumhist)]/length(logs_cumhist),nois(j) );
        else
            % AQUI DA PARA MELHORAR, ELE CAI AQUI QUANDO NAO H� ESTATISTICA  SUFICIENTE
            variable_uniform(i,j) = 0.5;
            logs_sub2 = nois(j);
        end
        
        index = abs(reference_variables_filtered(:,j) - logs_sub2) < grid_size; %filtra apenas valores em torno do valor amostrado
        reference_variables_filtered(~index,:) = NaN;
    end
end
% toc


variable_uniform(isnan(variable_uniform)) = 1e-5;
variable_uniform(variable_uniform==0) = 1e-5;
variable_uniform(variable_uniform==1) = 1-1e-5;


end