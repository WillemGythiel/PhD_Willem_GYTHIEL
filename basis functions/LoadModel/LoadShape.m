function [P] = LoadShape(Nodes,Trigs,Load,varargin)
if nargin<4
    Lfunc = @(Cor) Load_Shape3D(Cor);
    Ld_Nd = Load_Builder(Nodes,Trigs,Lfunc);
else
    dNodesdx = varargin{1};
    Lfunc = @(Cor) Load_Shape3D(Cor,'derivative');
    dNd = reshape(dNodesdx(:,2:4,:),3*size(Nodes,1),size(dNodesdx,3));
    Ld_Nd = Load_Builder(Nodes,Trigs,Lfunc,'derivative')*dNd;
end

nLc = size(Load,2);
P = zeros([size(Ld_Nd),nLc]);
for lc = 1:nLc
    Ld_Mag = diag(Load(:,lc));
    P(:,:,lc) = Ld_Mag*Ld_Nd;
end
end



function [P] = Load_Builder(Nodes,Trigs,Lfunc,varargin)
nNd = size(Nodes,1);    nLds = 6*nNd;
if nargin < 4
    nC = 1;
elseif strcmp(varargin{1},'derivative')
    nC = 3*nNd;
else
    print('Input not supported')
end
Pv = zeros(nLds,nC);
for i = 1:size(Trigs,1)
    NdT  =  Trigs(i,2:4);
    Cor  = Nodes(NdT,2:4);
    Li   = Addkron(6*(NdT-1),(1:6));
    if nargin < 4
        Lj = 1;
    else
        Lj = Addkron(NdT,nNd*(0:2));        
    end
    Lpn = Lfunc(Cor);
    Pv(Li,Lj) = Pv(Li,Lj)+Lpn;
end
[Pi,Pj] = ndgrid((1:nLds),(1:nC));
P = sparse(Pi(:),Pj(:),Pv(:),nLds,nC);
end




function [A] = Load_Shape3D(Cor,varargin)
x1 = Cor(1); y1 = Cor(4); z1 = Cor(7);
x2 = Cor(2); y2 = Cor(5); z2 = Cor(8);
x3 = Cor(3); y3 = Cor(6); z3 = Cor(9);
Ax = ((y1-y2)*(z1-z3)-(y1-y3)*(z1-z2));
Ay = ((z1-z2)*(x1-x3)-(z1-z3)*(x1-x2));
Az = ((x1-x2)*(y1-y3)-(x1-x3)*(y1-y2));
Norm = (Ax^2+Ay^2+Az^2)^(0.5);
if nargin<2
    A = Norm;
elseif strcmp(varargin{1},'derivative')
    dNorm = [((y2-y3)*Az-(z2-z3)*Ay)/Norm,(-(x2-x3)*Az+(z2-z3)*Ax)/Norm,((x2-x3)*Ay-(y2-y3)*Ax)/Norm,...
            (-(y1-y3)*Az+(z1-z3)*Ay)/Norm,((x1-x3)*Az-(z1-z3)*Ax)/Norm,(-(x1-x3)*Ay+(y1-y3)*Ax)/Norm,...
            ((y1-y2)*Az-(z1-z2)*Ay) /Norm,(-(x1-x2)*Az+(z1-z2)*Ax)/Norm,((x1-x2)*Ay-(y1-y2)*Ax)/Norm];
    A = dNorm;
else
    print('Input not supported')
end
A = repmat(A/6,18,1);
end

% function [A] = Load_ShapeProj(Cor,varargin)
% x1 = Cor(1); y1 = Cor(4); z1 = Cor(7);
% x2 = Cor(2); y2 = Cor(5); z2 = Cor(8);
% x3 = Cor(3); y3 = Cor(6); z3 = Cor(9);
% if nargin<2
%     Ax = ((y1-y2)*(z1-z3)-(y1-y3)*(z1-z2));
%     Ay = ((z1-z2)*(x1-x3)-(z1-z3)*(x1-x2));
%     Az = ((x1-x2)*(y1-y3)-(x1-x3)*(y1-y2));
% elseif strcmp(varargin{1},'derivative')
%     Ax = [0,z2-z3,-y2+y3,0,-z1+z3,y1-y3,0,z1-z2,-y1+y2];
%     Ay = [-z2+z3,0,x2-x3,z1-z3,0,-x1+x3,-z1+z2,0,x1-x2];
%     Az = [y2-y3,-x2+x3,0,-y1+y3,x1-x3,0,y1-y2,-x1+x2,0];
% else
%     print('Input not supported')
% end
% A = repmat([Ax;Ay;Az]/6,6,1);
% end

