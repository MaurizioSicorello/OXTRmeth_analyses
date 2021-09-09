cd('.\Data')

% load data
DF = readmatrix('MediationData.csv');
DF = zscore(DF);

xx = num2cell(DF(:,1));
yy = num2cell(DF(:,2));
xx_cat = num2cell(DF(:,3));
mm = DF(:, 4:size(DF, 2));
mm_cell = cell(size(xx,1), 1);

for i=1:size(DF, 1)
   
    mm_cell{i} = mm(i,:)';
    
end

%compute PDM for continuous CTQ
pdm = multivariateMediation(xx,yy,mm_cell, 'B', 20, 'svd', 'plots');

t = multivariateMediation(pdm, 'nPDM', 4, 'svd')

varRetained = pdm.dat.M_tilde*pdm.dat.Dt;
sum(var(varRetained)/size(mm, 2))

cd('..\Results')
pdmBoot = multivariateMediation(pdm, 'nPDM', 4, 'svd', 'bootPDM', 'Bsamp', 5000, 'save2file', 'multimedBootResults_cont.mat', 'returnbootsamples');


% compute PDM for dichotomous trauma
pdm2 = multivariateMediation(xx_cat,yy,mm_cell, 'B', 20, 'svd', 'plots');

varRetained = pdm2.dat.M_tilde*pdm2.dat.Dt;
sum(var(varRetained)/size(mm, 2))

% cd('..\Results')
% pdm2Boot = multivariateMediation(pdm, 'nPDM', 4, 'svd', 'bootPDM', 'Bsamp', 5000, 'save2file', 'multimedBootResults_cat.mat', 'returnbootsamples');


cd('Results')

modelCont = load('multimedBootResults_cont.mat')
bootABcont = modelCont.out.boot.SamplesTheta{1}(3,:).*modelCont.out.boot.SamplesTheta{1}(4,:);
modelCont.out.Theta{1}(4)
quantile(bootABcont, 0.025)
quantile(bootABcont, 0.975)
plot(modelCont.out.Wfull{1})

modelCat = load('multimedBootResults_cat.mat')
bootABcat = modelCat.out.boot.SamplesTheta{1}(3,:).*modelCat.out.boot.SamplesTheta{1}(4,:);
modelCat.out.Theta{1}(4)
quantile(bootABcat, 0.025)
quantile(bootABcat, 0.975)
plot(modelCat.out.Wfull{1})

t = table(modelCont.out.Wfull{1}, modelCat.out.Wfull{1},'VariableNames', {'Wcont', 'Wcat'})


length(modelCont.out.Wfull{1})

