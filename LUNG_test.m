clear;
clc;

addpath('./tool');
[biodata,surv]=readbiodata('LUNG');
biodata=dataproc(biodata);

opts.localK1=50;
opts.localK2=50;
opts.clusternum=4;
opts.alpha =1;    
[S,w] = AWSL(biodata,opts);

group = SpectralClustering(S,opts.clusternum);
