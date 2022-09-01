function [Perimeter] = get_perimeter(Nodes,Elements,domain)
ElemNd    =  unique(Elements(:,7));
vrts      = domain.Vertices;
nVrt      = size(vrts,1);
VertNr    = zeros(nVrt,1);
eq        = @(n1,n2) abs(n1(:,1)-n2(:,1))+abs(n1(:,2)-n2(:,2))<1e-3;
for i=1:nVrt
    VertNr(i) = Nodes(eq(Nodes(:,2:3),vrts(i,:)),1);
end
perm = [VertNr, circshift(VertNr,-1)];
Perimeter = cell(nVrt,2);
onLine    = @(n1,n2,nd) abs((n1(1)-n2(1))*(n1(2)-nd(:,2))-(n1(1)-nd(:,1))*(n1(2)-n2(2)))<1e-3;
for i=1:nVrt
    Perimeter{i,1} = perm(i,:);
    Perimeter{i,2} = setdiff(Nodes(onLine(Nodes(perm(i,1),2:3),Nodes(perm(i,2),2:3),Nodes(:,2:3)),1),[perm(i,:),ElemNd]);
end
end