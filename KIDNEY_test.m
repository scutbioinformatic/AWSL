clear;
clc;

addpath('./tool');
[biodata,surv]=readbiodata('KIDNEY');
biodata=dataproc(biodata);

opts.localK1=55;
opts.localK2=55;
opts.clusternum=3;
opts.alpha =1;    
[S,w] = AWSL(biodata,opts);

group = SpectralClustering(S,opts.clusternum);