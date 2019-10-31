function [S,w] = AWSL(X,opts)
%Auto-weighted similarity learning for Multi-omic data
%Input:
% X    -- cell, multi-omic data
% X{i} -- i-th matrix, row for feature, column for sample d*n
%opts.localK1 -- the neighbor number of z_i
%opts.localK2 --the neighbor number of x_i
%opts.clusternum  --cluster number
%opts.alpha  --control the  block diagonal matrix induced regularizer

%Output:
% S    --S Affinity Matrix n*n
% w -- weight of each omic data types.

num_omics = length(X);
n = size(X{1},2);

% setting default parameters
beta = 1;
alpha = 1;
clusternum=3;

%initial S S1 S2,and w;
S=zeros(n,n);
S1=S;
S2=S;
w=1/num_omics*ones(num_omics,1);

%initial parameters related to ADMM
num_iter = 200;
mu1=1e-2;
mu2=1e-2;
Lambda1 = zeros(n,n);
Lambda2 = zeros(n,n);
err_thr = 1e-5;
rho = 1.2;
max_mu = 1e6;


if ~exist('opts', 'var')
    opts = [];
else
    if ~isstruct(opts)
        error('Parameter error: opts is not a structure.');
    end
end

if isfield(opts, 'alpha');         alpha = opts.alpha;	end
if isfield(opts, 'clusternum');	   clusternum = opts.clusternum;end
if isfield(opts, 'localK1');	   localK1 = opts.localK1;end
if isfield(opts, 'localK2');	   localK2 = opts.localK2;end

iter = 0;

Am=ConstructAm(X,localK1,localK2);

while(iter < num_iter)
A=zeros(n,n);
for v=1:num_omics
    A=w(v)*Am{v}+A;
end
    iter=iter+1;
    %updata W
    D=diag(S2*ones(n,1));
    [U,~,~]=eig1(D-S2,clusternum,0);
    T=U*U';

    %updata S
    S=(mu1*(S1-Lambda1/mu1)+mu2*(S2-Lambda2/mu2)+A)/(mu1+mu2);
    
    %update S1
    for i=1:n
        S1(:,i)=EProjSimplex_new(S(:,i)+Lambda1(:,i)/mu1);
    end
    
    %update S2
    S2tmp=S+Lambda2/mu2+alpha/mu2*(diag(T)*ones(1,n)-T)-diag(diag(S+Lambda2/mu2+alpha/mu2*(diag(T)*ones(1,n)-T)));
    S2=max((S2tmp+S2tmp')/2,0);
    
    %update w
    zsz=zeros(num_omics,1);
    for v=1:num_omics
        zsz(v)=trace(Am{v}*S);
    end
    w=Hbeta(zsz,1/beta);
    
    %update multipliers and penalty scalars
    Lambda1=Lambda1+mu1*(S-S1);
    Lambda2=Lambda2+mu2*(S-S2);
    mu1=max(rho*mu1,max_mu);
    mu2=max(rho*mu2,max_mu);
    if(iter>=2 && max([norm(S-S1,'inf'),norm(S-S2,'inf')])<=err_thr)
        break
    end
end
end
