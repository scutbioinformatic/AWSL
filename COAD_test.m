clear;
clc;

addpath('./tool');
[biodata,surv]=readbiodata('COAD');
biodata=dataproc(biodata);

opts.localK1=90;
opts.localK2=90;
opts.clusternum=3;
opts.alpha =0.1;    
[S,w] = AWSL(biodata,opts);

group = SpectralClustering(S,opts.clusternum);