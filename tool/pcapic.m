function pcapic(S,group)
[COEFF SCORE latent]=pca(S);
pcaData1=SCORE(:,1:2);
%silvalue=mean(silhouette(ConZ,group));
%fprintf('%f',silvalue);
pcaData2=SCORE(:,1:3);
figure,scatter3(pcaData2(:,1), pcaData2(:,2), pcaData2(:,3), 10, group, 'filled');
figure,SIMLR_DisplayVisualization(pcaData1,group,215,[],10,800,0.1,'true')
end