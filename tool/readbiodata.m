function [biodata,surv]=readbiodata(type)
% sur: survial information
%biodata  row for sample, column for feature 
% type £ºKIDNEY,GBM,COLON,BIC, LUNG

switch type
    case 'GBM'
        importdata('.\biodata\GBM\GLIO_Gene_Expression.txt');
        X1=ans.data;
        importdata('.\biodata\GBM\GLIO_Methy_Expression.txt');
        X2=ans.data;
        importdata('.\biodata\GBM\GLIO_Mirna_Expression.txt');
        X3=ans.data;
        importdata('.\biodata\GBM\GLIO_survival.txt');
        surv=ans.data;
        biodata={X1,X2,X3};
    case 'BIC'
        importdata('.\biodata\BREAST\BREAST_Gene_Expression.txt');
        X1=ans.data;
        importdata('.\biodata\BREAST\BREAST_Methy_Expression.txt');
        X2=ans.data;
        importdata('.\biodata\BREAST\BREAST_Mirna_Expression.txt');
        X3=ans.data;
        importdata('.\biodata\BREAST\BREAST_survival.txt');
        surv=ans.data;
        biodata={X1,X2,X3};
end
end
