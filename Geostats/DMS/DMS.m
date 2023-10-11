function [logs_simulated_all] = DMS(I,J, range, type, angles, grid_size, ref_variables, cond_pos, cond_value, num_of_sims)
% TO DO: EXPLAIN INPUTS


krig_mean = zeros(size(ref_variables,2), I, J);
krig_std = ones(size(ref_variables,2), I, J);
%% if we have conditionin points, use kriging to condition FFTMA simulation
if size(cond_value, 1) > 0

    % Gaussian transforming the conditional points
    cond_value_uniform = nonParametric_to_uniform( cond_value, ref_variables , grid_size);        
    cond_value_gaussian_2krig = norminv(cond_value_uniform);
    
    [X,Y] = meshgrid(1:I,1:J);
    xcoords = [ Y(:) X(:)];
    for i = 1:1:size(krig_mean,1)
        krig = 1;
        [mean_krig, var_krig] = Kriging_options(xcoords, cond_pos, cond_value_gaussian_2krig(:,i), 1, range, type, krig, angles);
        krig_mean(i,:) = mean_krig(:);
        krig_std(i,:) = sqrt(var_krig(:));
    end
end

% correlation function to use FFTMA in DMS. It is also possible to use SGS instead
[correlation_function] = construct_correlation_function(zeros(I,J), range, type, angles);
       
%% DMS 
logs_simulated_all = cell(num_of_sims,1);
parfor n = 1:1:num_of_sims
    
    % simulating Gaussian realizations using FFTMA
    simulations2D = zeros(size(ref_variables,2), I, J);
    for i = 1:1:size(ref_variables, 2)
        white_noise = randn(I,J);
        simulations_gaussian = FFT_MA_3D(correlation_function,white_noise);
        simulations_gaussian = make_it_gaussian(simulations_gaussian(:)); % It makes the simulation to be "Perfectly" Gaussian distributed, therefore, better non-parametric simulations
        simulations_gaussian = reshape(simulations_gaussian,I,J);
        simulations2D(i,:,:) = simulations_gaussian;
    end
    
    simulations_gaussian = simulations2D(:,:)';
    simulations_conditioned = zeros(size(simulations_gaussian));
    
    % conditioning Gaussian realizations using kriging results
    for i = 1:1:size(simulations2D,1)
        mean_krig = reshape(krig_mean(i,:), size(simulations2D,2), size(simulations2D,3));
        var_krig = reshape(krig_std(i,:), size(simulations2D,2), size(simulations2D,3));
        simulation_gaussian = reshape(simulations_gaussian(:,i), size(simulations2D,2), size(simulations2D,3));
        simulation_gaussian = simulation_gaussian.*var_krig +  mean_krig;
        simulations_conditioned(:,i) = reshape(simulation_gaussian, I * J, 1);
    end
    
    % Transforming the conticioned Gaussian realizations to the non parametris pdf
    simulations_conditioned = normcdf(simulations_conditioned);
    simulations_conditioned_nonParametric = uniform_to_nonParametric( simulations_conditioned, ref_variables, grid_size);
    
    result = zeros(size(simulations_conditioned_nonParametric,2), I, J);
    for i = 1:1:size(simulations_conditioned_nonParametric,2)
        result(i,:,:) = reshape(simulations_conditioned_nonParametric(:,i), I, J);
    end
    logs_simulated_all{n} = result;
end


