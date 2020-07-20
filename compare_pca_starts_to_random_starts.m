% This script shows an instance where using random restarts for SigClust gives a
% significant result, while PCA starts do not.

%% Read TCGA kidney cancer data (n=50)
kid = readmatrix('kidney_cancer_data.csv');

%% Run SigClust with PCA starts
rng(123)
paramstruct = struct('vclass',0, ...
                     'iCovEst',1, ...
                     'nsim',500, ...
                     'iBGSDdiagplot',0, ...
                     'iCovEdiagplot',0, ...
                     'ipValplot',1, ...
                     'iscreenwrite',2, ...
                     'twoMtype', 2) ;  % PCA starts 
[pval_pca_starts, zscore_pca_starts] = SigClustSM(kid, paramstruct);
saveas(gcf, 'pca_starts.png');
%% Random starts
rng(123)
clf;
paramstruct = struct('vclass',0, ...
                     'iCovEst',1, ...
                     'nsim',500, ...
                     'iBGSDdiagplot',0, ...
                     'iCovEdiagplot',0, ...
                     'ipValplot',1, ...
                     'iscreenwrite',2, ...
                     'twoMtype', 1) ;  % Random starts 
[pval_random_starts, zscore_random_starts] = SigClustSM(kid, paramstruct);
saveas(gcf, 'random_starts.png');

%% PCA starts vs random starts for random matricies
% Warning: I ran this overnight!
max_exponent = 14;
rng(123)
results = zeros(max_exponent, 2*max_exponent);
for i = (1:max_exponent)
    disp({'i', i})
    for j = (1:max_exponent)
        disp({'j', j})
        grm = randn(2^i, 2^j);
      
        [bestclass, vindex, midx] = SigClust2meanRepSM(grm);
        results(i, 2*j-1) = min(vindex);
        
        [bestclass, bestCI] = SigClust2meanFastSM(grm, struct('ioutplot', 0));
        results(i, 2*j) = bestCI;
    end
end

save('results.mat', 'results');

%% Create heatmap
ylabs = ["2^1", "2^2", "2^3", "2^4", "2^5", "2^6", "2^7", "2^8", "2^9", "2^{10}", "2^{11}", "2^{12}", "2^{13}", "2^{14}"];

xlabs = reshape([ylabs+"R"; ylabs+"P"], size(ylabs,1), []);

heatmap(xlabs, ylabs, results, ...
       'XLabel', 'Replicates (R = Random, P = PCA)', ...
       'YLabel', 'Dimension', ...
       'Title', '2 means clustering performance on N(0,1) matrices, random vs PCA');

%% We confirm that the simulated CIs for both methods are the same,
% meaning that the issue is in the computation of the data cluster index

kid_T = kid';
d = size(kid_T,1) ;
n = size(kid_T,2) ;
outstruct = pcaSM(kid_T) ;
veigval = getfield(outstruct,'veigval') ;
nev = length(veigval) ;
veigval = [veigval; zeros(d-nev,1)] ;
vdata = reshape(kid_T, d*n, 1) ;
simbackvar = madSM(vdata)^2;
vsimeigval = SigClustCovEstHH(veigval, simbackvar) ;
vscale = sqrt(vsimeigval) ;
vSimIndex_Random = [] ;
vSimIndex_PCA = [];
twoMsteps = 1;

for isim = 1:500 
    mdatsim = randn(d,n) ;
    mdatsim = mdatsim .* vec2matSM(vscale,n) ;
    
    % Compute CI with random (k-means++) start
    paramstruct = struct('nrep',twoMsteps) ;
    [~, vindex] = SigClust2meanRepSM(mdatsim,paramstruct) ;
    ClustIndSim_Random = min(vindex) ; % SigClust2meanRepSM returns full vector of Clust Indices

    % Compute CI with PCA start
    paramstruct = struct('maxstep',twoMsteps,  ...
                         'ioutplot',0) ;
    [~, ClustIndSim_PCA]= SigClust2meanFastSM(mdatsim,paramstruct) ;

    vSimIndex_Random = [vSimIndex_Random; ClustIndSim_Random] ;
    vSimIndex_PCA = [vSimIndex_PCA; ClustIndSim_PCA] ;
end
    
scatter(vSimIndex_Random, vSimIndex_PCA)
xlabel('Random Start')
ylabel('PCA Start')
title("Null cluster indices, PCA vs random starts")
saveas(gcf, 'scatter_pca_vs_random.png')

