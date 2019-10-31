function alldata=dataproc(biodata)
Geneexpression = biodata{1};
Geneexpression(Geneexpression>10) = 10;
Geneexpression(Geneexpression<-10) = -10;

Methyexpression = biodata{2};

rnaexpression = biodata{3};
rnaexpression(rnaexpression>10) = 10;
rnaexpression(rnaexpression<-10) = -10;

Data1 = Standard_Normalization(Geneexpression');
Data2 = Standard_Normalization(Methyexpression');
Data3 = Standard_Normalization(rnaexpression');
alldata={Data1',Data2',Data3'};
end