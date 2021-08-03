% load data
DF = readmatrix('MediationData.csv');
DF = zscore(DF);

xx = num2cell(DF(:,1));
yy = num2cell(DF(:,2));
mm = DF(:, 3:size(DF, 2));
mm_cell = cell(size(xx,1), 1);

for i=1:size(DF, 1)
   
    mm_cell{i} = mm(i,:)';
    
end

% compute PDMs
pdm = multivariateMediation(xx,yy,mm_cell, 'B', 50, 'svd', 'plots');

% bootstrap first four PDMs
% pdm2 = multivariateMediation(pdm, 'nPDM', 4, 'svd', 'bootPDM', 'Bsamp', 5000, 'save2file', 'multimedBootResults.mat', 'returnbootsamples');

% load results
bootResults = load('multimedBootResults.mat');

% plot bootstrapping estimates of first theta, last entry
histogram(bootResults.out.boot.SamplesTheta{1}(4,:))
