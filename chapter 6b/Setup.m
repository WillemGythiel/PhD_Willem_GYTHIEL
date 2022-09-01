function Problem = Setup(Model)
Problem = struct;
Problem.R = 10;     % Radius of length scale
Problem.int =  get_intNodes(Model.Nodes,Model.Elements,Model.info);
gradient = 'both';
%--------------------------------------------------------shape optimization
if any(strcmp(gradient,'shape'))
   Problem.BID = []; xminS = []; xmaxS = []; xvalS = [];
   [Problem,xvalC,xminC,xmaxC] = setup_shape(Model,Problem);
elseif any(strcmp(gradient,'size'))
    [Problem,xvalS,xminS,xmaxS] = setup_size(Model,Problem);
    xminC = []; xmaxC = []; xvalC = [];
    Problem.dNodesdx = zeros([size(Model.Nodes),size(Problem.BID,1)]);
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
int            = get_intNodes(Model.NodesCon,Model.Elements,Model.info);
Problem.CID    = Addkron(int(:,1),[1;2;3]*size(Model.NodesCon,1));
Problem.SID    = Addkron(3*(int(:,1)-1),[1;2;3]);
Problem.nCvar  = size(Problem.CID,1);
Mini = @(dim) min(Model.NodesCon(:,dim)); 
Maxi = @(dim) max(Model.NodesCon(:,dim)); 
xminC = Addkron(zeros(size(int,1),1),[Mini(2);Mini(3);Mini(4)]);
xmaxC = Addkron(zeros(size(int,1),1),[Maxi(2);Maxi(3);Maxi(4)*2]);
xvalC = Model.NodesCon(Problem.CID);
step = Model.info.stepP;
for desVar = 1:Problem.nCvar
    if mod(Problem.SID(desVar),3)==0
        xminC(desVar) = max(xminC(desVar),xvalC(desVar)-5);
        xmaxC(desVar) = min(xmaxC(desVar),xvalC(desVar)+5);
    else
        xminC(desVar) = max(xminC(desVar),xvalC(desVar)-step);
        xmaxC(desVar) = min(xmaxC(desVar),xvalC(desVar)+step);
    end
end
end

function [Problem,xvalS,xminS,xmaxS] = setup_size(Model,Problem)
Elem = Model.Elements;
Problem.BID   = unique(Elem(:,3));
nS    = size(Problem.BID,1);
xminS = ones(nS,1)*1e-3;
xmaxS = ones(nS,1)+1e-10;
xallS = Model.Sections(:,8)*2/Model.cache.h; 
xvalS = xallS(Problem.BID);
end

function [NdInt] = get_intNodes(Nodes,Elements,info)
Perimeter = get_perimeter(Nodes,Elements,info.domain);
NdVrt = unique(cell2mat(Perimeter(:,1)));
NdEdg = cell2mat(Perimeter(:,2));
NdInt = Nodes(setdiff(Nodes(:,1),[unique(Elements(:,7));NdVrt;NdEdg]),:);
end
