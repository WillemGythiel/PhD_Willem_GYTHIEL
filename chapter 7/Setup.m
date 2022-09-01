function [Problem] = Setup(Model)
Problem = struct;
Problem.R = 10;     % Radius of length scale
err = 1e-2;
gradient = 'both';
%--------------------------------------------------------shape optimization
Problem.CID = []; Problem.SID = []; xminC = []; xmaxC = []; xvalC = [];
Node = Model.Node;
num   = size(Node,1);
NdID = Node(1:end-1,1);
Per = @(d) abs(Node(:,d)-17/2)>17/2-err;
Dia = NdID(~Per(2)&~Per(3)& abs(Node(:,2)-Node(:,3))<err);
int = NdID(~Per(2)&~Per(3)& abs(Node(:,2)-Node(:,3))>=err);
Problem.DID = Dia;  % Diagonal nodes
Problem.int = int;
if any(strcmp(gradient,{'shape','both'}))
    PrX = NdID(Per(2)&~Per(3));
    PrY = NdID(Per(3)&~Per(2));
    Problem.CID = [Addkron(int,[1;2;3]*num);PrY+1*num;PrY+3*num;PrX+2*num;Dia+1*num;Dia+3*num];  % Indices of xVars in Nodes
    Problem.SID = [Addkron(3*(int-1),[1;2;3]);3*(PrY-1)+1;3*(PrY-1)+3;3*(PrX-1)+2;3*(Dia-1)+1;3*(Dia-1)+3]; % Indices of xVars in xAll list
    nC    = size(Problem.CID,1);
    xminC = -17*ones(nC,1);
    xmaxC =  17*ones(nC,1);
    xminC(mod(Problem.SID,3)==0)=0;
    xmaxC(mod(Problem.SID,3)==0)=4.0;
    xvalC = Node(Problem.CID);
end
%---------------------------------------------------------size optimization
Problem.BID = []; xminS = []; xmaxS = []; xvalS = [];
if any(strcmp(gradient,{'size','both'}))
    Elem = Model.Elem;
    Problem.BID   = unique(Elem(:,3));
    nS    = size(Problem.BID,1);
    xminS = ones(nS,1)*1e-3;
    xmaxS = ones(nS,1)+1e-2;
    xallS = Model.Sections(:,8)*2/Model.cache.h; 
    xvalS = xallS(Problem.BID);
end
%--------------------------------------------------------------------values
Problem.lb = [xminC;xminS];               
Problem.ub = [xmaxC;xmaxS];
Problem.x  = [xvalC;xvalS];
Problem.SigmaT = 355e3;
Problem.SigmaC = 355e3;
Problem.nCvar = size(xminC,1);
Problem.factor = 600;
end
