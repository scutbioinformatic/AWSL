function Am=ConstructAm(X,localK1,localK2)
num_omics = length(X);
n = size(X{1},2);
localK1=min(localK1,n);
alpha=0.5;
for v=1:num_omics
    dist=dist2(X{v}',X{v}');
    Wall = affinityMatrix(dist, localK1, alpha);
    Wall = Wall./repmat(sum(Wall,2),1,n);
    Wall = (Wall + Wall')/2;
    Z{v}=FindDominateSet(Wall,localK1);
end
localK2=min(localK2,n-1);
for v=1:num_omics
    distance=dist2(Z{v}',Z{v}');
    [sorted,index] = sort(distance);
    neighborhood{v}= index(2:(1+localK2),:);
end
Am=repmat({zeros(n,n)},1,num_omics);
for v=1:num_omics
    for i=1:n
        zneigh=repmat(Z{v}(:,i),1,localK2)-Z{v}(:,neighborhood{v}(:,i));
        Am{v}=zneigh*zneigh'+Am{v};
    end
end
end
function newW = FindDominateSet(W,K)
[m,n]=size(W);
[YW,IW1] = sort(W,2,'descend');
clear YW;
newW=zeros(m,n);
temp=repmat((1:n)',1,K);
I1=(IW1(:,1:K)-1)*m+temp;
newW(I1(:))=W(I1(:));
newW=newW./repmat(sum(newW,2),1,n);
clear IW1;
clear IW2;
clear temp;
end

