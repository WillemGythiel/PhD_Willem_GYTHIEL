function Problem = Setup(Model)
Problem = struct;
Problem.R = 10;     % Radius of length scale
Problem.int =  get_intNodes(Model.Node(1:end-1,:),Model.Elements,Model.info);
gradient = 'size';
%--------------------------------------------------------shape optimization
if any(strcmp(gradient,'shape'))
   Problem.BID = []; xminS = []; xmaxS = []; xvalS = [];
   [Problem,xvalC,xminC,xmaxC] = setup_shape(Model,Problem);
elseif any(strcmp(gradient,'size'))
    [Problem,xvalS,xminS,xmaxS] = setup_size(Model,Problem);
    xminC = []; xmaxC = []; xvalC = [];
    Problem.dNode    = zeros([size(Model.Node),size(Problem.BID,1)]);
    Problem.NodL     = [];
    Problem.nCvar    = 0;
elseif any(strcmp(gradient,'both'))
    [Problem,xvalS,xminS,xmaxS] = setup_size(Model,Problem);
    [Problem,xvalC,xminC,xmaxC] = setup_shape(Model,Problem);
end
%--------------------------------------------------------------------values
Problem.lb = [xminC;xminS];               
Problem.ub = [xmaxC;xmaxS];
Problem.x  = [xvalC;xvalS];
Problem.SigmaT = 350e3;
Problem.SigmaC = 350e3;
Problem.getPeri = @(Nodes,Elements,domain) get_perimeter(Nodes,Elements,domain);
end

function [Problem,xvalC,xminC,xmaxC] = setup_shape(Model,Problem)
[dNode,NodL,t,nCvar] = Make_dNodesdx(Model.Node,Model.Elements,Model.info,Problem.BID);
Problem.dNode    = dNode;
Problem.NodL     = NodL;
Problem.nCvar    = nCvar;
xvalC = max(t,1e-10);
xminC = zeros(nCvar,1);
xmaxC = ones(nCvar,1)+1e-10;    
end


function [Problem,xvalS,xminS,xmaxS] = setup_size(Model,Problem)
Elem = Model.Elements;
Problem.BID   = unique(Elem(:,3));
nS    = size(Problem.BID,1);
xminS = ones(nS,1)*1e-2;
xmaxS = ones(nS,1)+1e-10;
xallS = Model.Sections(:,8)/Model.cache.R; 
xvalS = xallS(Problem.BID);
end

function [NdInt] = get_intNodes(Nodes,Elements,info)
Perimeter = get_perimeter(Nodes,Elements,info.domain);
NdVrt = unique(cell2mat(Perimeter(:,1)));
NdEdg = cell2mat(Perimeter(:,2));
NdInt = Nodes(setdiff(Nodes(:,1),[unique(Elements(:,7));NdVrt;NdEdg]),:);
end


function [dNodesdx,NodL,t,nCvar] = Make_dNodesdx(Nodes,Elements,info,BID)
suppXY = cell2mat(info.suppNode(:,1));
vertXY = info.domain.Vertices;
check  = [suppXY(:,1:2);vertXY];
nCons  = size(check,1);
remove = zeros(nCons,1);
eq     = @(n1,n2) abs(n1(:,1)-n2(:,1))+abs(n1(:,2)-n2(:,2))<1e-5;
for i=1:nCons
    remove(i) = min(Nodes(eq(Nodes(:,2:3),check(i,:)),1));
end
remove2 = Nodes(abs(Nodes(:,4))>1e3,1);
Perimeter = get_perimeter(Nodes(1:end-1,:),Elements,info.domain);
% matrices with special nodes
NdVrt = unique(cell2mat(Perimeter(:,1)));
NdEdg = cell2mat(Perimeter(:,2));
NdInt = Nodes(setdiff(Nodes(:,1),[unique(Elements(:,7));remove;remove2;NdEdg]),:);
% just some numbers
nNdEdg  = numel(nonzeros(cell2mat(Perimeter(:,2))));
nCvar   = 2*nNdEdg+3*size(NdInt,1);
nDesVar = nCvar+size(BID,1);
% making dNodesdx
dNodesdx = zeros([size(Nodes),nDesVar]);

t = zeros(nCvar,1);
NodL = Nodes;

NdAux = [NdVrt;NdEdg;NdInt(:,1)];
dvAux = @(dim) max(Nodes(NdAux,dim))-min(Nodes(NdAux,dim));
dv    = [0,dvAux(2),dvAux(3),max(dvAux(4),1e-10)];
tv = @(Ni,dim) (Ni(dim)-min(Nodes(NdAux,dim)))/dv(dim);
% dv    = [0,dvAux(2),dvAux(3),4];

dvP = @(nd,dim) Nodes(nd(2),dim)-Nodes(nd(1),dim);
tvP = @(nE,nV,dim) (Nodes(nE,dim)-Nodes(nV(1),dim))/(Nodes(nV(2),dim)-Nodes(nV(1),dim));
counter = 1;
for i=1:numel(NdVrt)
    NdExtrI = Perimeter{i,1};
    NdEdgeI = Perimeter{i,2};
    for j=1:numel(NdEdgeI)
        dNodesdx(NdEdgeI(j),2:3,counter) = [dvP(NdExtrI,2),dvP(NdExtrI,3)];
        dNodesdx(NdEdgeI(j),4,counter+1) = dv(4);
        NodL(NdEdgeI(j),2:4) = [Nodes(NdExtrI(1),2:3),min(Nodes(NdAux,4))];
        t(counter+(0:1)) = [min([tvP(NdEdgeI(j),NdExtrI,2),tvP(NdEdgeI(j),NdExtrI,3)]);tv(Nodes(NdEdgeI(j),:),4)];
        counter = counter+2;
    end
end

for i=1:size(NdInt,1)
    dNodesdx(NdInt(i),2:4,counter+(0:2)) = reshape(diag(dv(2:4)),1,3,3);
    NodL(NdInt(i),2:4) = [min(Nodes(NdAux,2)),min(Nodes(NdAux,3)),min(Nodes(NdAux,4))];
    t(counter+(0:2)) = [tv(NdInt(i,:),2),tv(NdInt(i,:),3),tv(NdInt(i,:),4)];
    counter = counter+3;
end
end