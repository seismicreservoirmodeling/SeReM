function [ noise3 ] = make_it_gaussian(noise)
%% Makes the noise "perfectly" gaussian distributed

%     noise = randn(500,1); %gera gaussiana empirica
    nn = noise;
    % figure;
    % hist(noise,50);
    MAX = max(noise);
    MIN = min(noise);
    noise = noise - min(noise);
    noise = noise ./ max(noise);


    %uniformiza
    logs_cumhist = unique(sort(noise));
    logs_cumhist_x = [1:1:size(logs_cumhist,1)]' -1; 
    logs_cumhist_x = logs_cumhist_x ./ max(logs_cumhist_x);
    noise2 = interp1(logs_cumhist, logs_cumhist_x, noise, 'pchip'); %amostra

    % figure;
    % hist(noise2, 50)

    x = (-4):0.001:4;
    gaussian = gaussmf(x,[1 0])';
    gaussian = unique(cumsum(gaussian));
    gaussian = gaussian - min(gaussian);
    gaussian = gaussian ./ max(gaussian);
%     figure
%     plot( gaussian);

    logs_cumhist = gaussian;
    logs_cumhist_x = [1:1:size(logs_cumhist,1)]' -1; 
    logs_cumhist_x = logs_cumhist_x ./ max(logs_cumhist_x);
    noise3 = interp1(logs_cumhist, logs_cumhist_x, noise2, 'pchip'); %amostra

    % figure; 
    % plot(nn)
    % hold on;
    % plot(noise3)

    noise3 = noise3 - mean(noise3);
    noise3 = noise3 ./ std(noise3);
end
% figure;
% hist(noise3, 50)
% 
% figure;
% plot(sort(noise3))
