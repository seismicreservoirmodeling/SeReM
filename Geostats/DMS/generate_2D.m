function [] = generate_2D(models,cond_pos)

figure

for var = 1:size(models,1)
    subplot(2,round(size(models,1)/2),var)
    imagesc(squeeze(models(var,:,:)));
    title('z^'+string(var));
   
    if nargin>1
        hold on
        scatter(cond_pos(:,2),cond_pos(:,1), 'black','x');
    end
    
    
end