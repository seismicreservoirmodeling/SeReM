function [reference_variables] = extend_dateset_KDE(reference_variables,n_times,std_normalized)

min2norm = min(reference_variables);

reference_variables = reference_variables - repmat(min2norm, size(reference_variables,1), 1);
max2norm = max(reference_variables);
reference_variables = reference_variables ./ repmat(max2norm, size(reference_variables,1), 1);

reference_variables_ = [];
for sampling = 1:n_times-1
    reference_variables_ = [reference_variables_ ; reference_variables + std_normalized*randn(size(reference_variables))];
end
reference_variables = [reference_variables ; reference_variables_];

reference_variables = reference_variables .* repmat(max2norm, size(reference_variables,1), 1);
reference_variables = reference_variables + repmat(min2norm, size(reference_variables,1), 1);