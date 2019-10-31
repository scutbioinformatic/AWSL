clear;
clc;

addpath('./tool');
[biodata,surv]=readbiodata('BIC');
biodata=dataproc(biodata);

opts.localK1=25;
opts.localK2=25;
opts.clusternum=5;
opts.alpha =0.01;    
[S,w] = AWSL(biodata,opts);

group = SpectralClustering(S,opts.clusternum);
