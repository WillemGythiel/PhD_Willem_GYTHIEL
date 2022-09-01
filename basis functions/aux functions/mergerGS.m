function [Node,Elem,Trig,xSec] = mergerGS...
         (Node,Elem,Trig,xSec,Dmin,Vertices,Perimeter)
[Node,ic] = mergeNodes(Node,Dmin,Vertices,Perimeter);
[Elem,xSec]  = mergeElements(Elem,xSec,ic);
[Trig]    = mergeTrigs(Trig,ic);
end
%% AUXILIAIRY FUNCTIONS %%
function [Node,ic] = mergeNodes(Node,Dmin,Vertices,Perimeter)
Z = linkage(Node(:,2:4),'complete');
% T = cluster(Z,'cutoff',Dmin,'criterion','distance');
T = cluster(Z,'MaxClust',size(Node,1)-1);
for i = 1:max(T)
    ids = Node(T==i,1);
    verts = intersect(ids,Vertices);
    perim = intersect(ids,Perimeter);
    if size(verts,1)==0 && size(perim,1)==0
        Node(ids,2:4) = repmat(mean(Node(ids,2:4),1),size(ids,1),1);
    elseif size(verts,1)==1
        Node(ids,2:4) = repmat(Node(verts,2:4),size(ids,1),1);
    elseif size(perim,1)>=1
        Node(ids,2:4) = repmat(mean(Node(perim,2:4),1),size(ids,1),1);
    else
        disp('Error: Distance between restricted nodes exceeds merge radius');
    end
end
[Ns,~,ic] = unique(round(Node(:,2:4),3),'rows','stable');
      nID = (1:size(Ns,1))';
    Node = [nID Ns];
end
function [Elem,xSec] = mergeElements(Elem,xSec,ic)
NdElem         = [ic(Elem(:,5)),ic(Elem(:,6)),ic(Elem(:,7))];
[NdElem,ja,jc]  = unique(sort(NdElem,2,'ascend'),'rows');
nelm        = size(NdElem,1);
xNew = zeros(nelm,1);
ElemNew        = Elem(ja,:);
ElemNew(:,5:7) = NdElem;
for i = 1:nelm
    xNew(i)             = sum(xSec(jc==i).^3)^(1/3);
end
keep = ElemNew(:,5)~=ElemNew(:,6);
Elem = ElemNew(keep,:);
xSec = xNew(keep);
Elem(:,1) = (1:size(Elem,1))';
Elem(:,3) = Elem(:,1);
end
function [Trig] = mergeTrigs(Trig,ic)
T1 = ic(Trig(:,2));
T2 = ic(Trig(:,3));
T3 = ic(Trig(:,4));
Trig(:,2:4) = [T1,T2,T3];
Trig        = Trig(T1~=T2 & T2~=T3 & T3~=T1,:);
Trig(:,1)   = (1:size(Trig,1))';
end
