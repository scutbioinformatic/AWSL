clear;
clc;

addpath('./tool');
[biodata,surv]=readbiodata('GBM');
biodata=dataproc(biodata);

opts.localK1=25;
opts.localK2=40;
opts.clusternum=3;
opts.alpha =10;    
[S,w] = AWSL(biodata,opts);

group = SpectralClustering(S,opts.clusternum);



