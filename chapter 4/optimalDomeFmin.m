function [Wgt] = optimalDomeFmin(radius)
load('DomeSubs.mat','domes');
Wgt = zeros(3,size(domes,1));
Wgt(1,:) = domes(:,1)';
for i=1:size(domes,1)
    if domes(i,1)==1
        [Nodes,Elements,dofs,P] = domeSchwedler(domes(i,2),domes(i,3),radius);
    elseif domes(i,1)==2
        [Nodes,Elements,dofs,P] = domeKiewitt(domes(i,2),domes(i,3),radius);
    else
        [Nodes,Elements,dofs,P] = domeFuller(domes(i,2),radius);
    end
   [weight,Ltot] = Optimizer(Nodes,Elements,dofs,P);
    Wgt(3,i) = weight;
    Wgt(2,i) = Ltot;
end
figure
hold on
scatter(Wgt(2,Wgt(1,:)==1),Wgt(3,Wgt(1,:)==1),'o');
scatter(Wgt(2,Wgt(1,:)==2),Wgt(3,Wgt(1,:)==2),'*');
scatter(Wgt(2,Wgt(1,:)==3),Wgt(3,Wgt(1,:)==3),'x');
hold off
xlabel('Total bar length') 
ylabel('Weight') 
end

function [weight,Ltot] = Optimizer(Nodes,Elements,dofs,P)
%% OPTIMIZATION %%
De   = 0.05;
Elm1 = ones(size(Elements(:,1)));
% Elm1 = 1;
x    = Elm1*De;
lb   = Elm1*0.0001;
ub   = Elm1*0.5;
x0   = x;
lts = sum((Nodes(Elements(:,6),2:4)-Nodes(Elements(:,5),2:4)).^2,2).^0.5;

fun     = @(x) Obj(Nodes,Elements,x);
nonlcon = @(x) Con(Nodes,Elements,dofs,P,x,lts);     
options = optimoptions('fmincon','StepTolerance',1e-10,'Algorithm','sqp',...
                  'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
                  'MaxIterations',50,'display','iter');

[~,fval,~,~] = fmincon(fun,x0,[],[],[],[],lb,ub,nonlcon,options);

weight = 100*fval*9.81;
Ltot   = sum(lts,1);

end


function [varargout] = Con(Nodes,Elements,dofs,P,x,lts)
%% FORMULATION OF THE OPTIMIZATION PROBLEM
fy = 235e6; sMax = fy/1.15; rho = 7850; Kf = 1; E = 210e9;
% Defining the initial, maximal and minamel design values
       Elm1 = ones(size(Elements(:,1)));
k    = Elm1*inf;
Types=  {1    'truss'};
Materials=  [1 E 0.3 rho];     
    % UPDATE SECTIONS MATRICES
    De = x;
    f  = 1-2*0.04;           
    I  = (De/2).^4*(1-f^4)*pi/2;
    A  = (De/2).^2*(1-f^2)*pi;
    Sections = [(1:size(x,1))' A k  k  I I/2 I/2];
%     Sections = [(1:size(x,1))' A inf  inf  I I/2 I/2];
    dI = 4*(De/2).^3*(1-f^4)*pi/2;
    dA = 2*(De/2)   *(1-f^2)*pi;
    Fmax    = (pi^2*E*I/2)./((Kf*lts).^2);
    dFmax = diag((pi^2*E*dI/2)./((Kf*lts).^2));
    fmax    = sMax*A; 
    dfmax = diag(sMax*dA);
    % NODE DEFINITION DERIVATIVES
    dNodes = zeros([size(Nodes),size(x,1)]);
    % SECTION DEFINITION DERIVATIVES
    dSections = zeros([size(Sections),size(x,1)]);
    dSections(:,2,:) = diag(dA);
    dSections(:,5,:) = diag(dI);
    dSections(:,6,:) = diag(dI)/2;
    dSections(:,7,:) = diag(dI)/2;
    % OUTPUT DISPLACEMENT AND SENSITIVITY
    [K,~,dK] = asmkm(Nodes,Elements,Types,Sections,Materials,dofs,dNodes,dSections);

    % Kiewitt dome
    U = K\P;
    FdF = zeros(size(dofs,1),size(x,1));
    for i = 1:size(x,1)
        FdF(:,i) = -dK{i}*U;
    end
    dU = zeros(size(U,1),1,size(x,1));
    dU(:,1,:) = K\FdF;

      Bars = [Elements(:,[1,5,6,7]),A,0*A+E];
      dof2 = round(3*(floor(dofs)-1)+100*mod(dofs,1),1);
      U2 = sparse(dof2,ones(size(dof2,1),1),U,3*size(Nodes,1),1);
      n2 = size(dU,3);
      dU2 = sparse(repmat(dof2,n2,1),Addkron((0:n2-1)',ones(size(dof2,1),1)),dU(:),3*size(Nodes,1),n2);
     N = Truss_F(Nodes,Bars,U2);
     dN = Truss_F(Nodes,Bars,U2,dU2,dA);
    f = [-N./Fmax-1;
          N./fmax-1];
    dfdx = [-(dN.*Fmax-dFmax.*N)./(Fmax.^2);
             (dN.*fmax-dfmax.*N)./(fmax.^2)]; 
    varargout{1}=f;
	varargout{2}=[]; varargout{3}=dfdx'; varargout{4}=[];
end

function [varargout] = Obj(Nodes,Elements,x)
%% FORMULATION OF THE OPTIMIZATION PROBLEM
rho = 7850; 
lts = sum((Nodes(Elements(:,6),2:4)-Nodes(Elements(:,5),2:4)).^2,2).^0.5;
    % UPDATE SECTIONS MATRICES
    De = x;
    f  = 1-2*0.04;           
    A  = (De/2).^2*(1-f^2)*pi;
    dA = 2*(De/2)   *(1-f^2)*pi;
    W       =  A'*lts*rho;
    dW      = dA.*lts*rho;
     f0   =  W/100;
    df0dx = dW/100;
    varargout{1}=f0;  varargout{2}=df0dx;

end


function [F] = Truss_F(Nodes,Bars,U,varargin)
nB  = size(Bars,1);
CsB = @(x) Nodes(Bars(x,2:4),2:4);
Ue  = @(x) U(Addkron(3*(Bars(x,2:3)-1),(1:3)));
F   =  zeros(nB,1);
Ke = @(x) KeL_Truss(CsB(x),Bars(x,5),Bars(x,6));
Te = @(x) Te_Truss(CsB(x));
if nargin < 4
    ind = @(x) 0*x+1;
    Us  = @(x) Ue(x);
    Qe  = @(x) 0;
else
    dU  = varargin{1};
    ind = @(x) x;
    Sec = varargin{2};
    Us  = @(x) dU(Addkron(3*(Bars(x,2:3)-1),(1:3)),:);
    dKe = @(x,y) KeL_Truss(CsB(x),Sec(x,:),Bars(x,6));
    Qe  = @(x) getIt(dKe(x)*Te(x)*Ue(x),2,1);
    F   = zeros(nB,nB);
end
for i = 1:nB
    Ftt = Ke(i)*Te(i)*Us(i);
    F(i,:)=Ftt(2,:);
    F(i,ind(i))=F(i,ind(i))+Qe(i);
end
end


%% Aux %%
function [el] = getIt(M,i,j)
el = M(i,j);
end

function [Te] = Te_Truss(Coor)
Xt = (Coor(2,:)-Coor(1,:));
Te  = kron(eye(2),Xt/norm(Xt));
end
function [KeL] = KeL_Truss(Coor,A,E)
L  = norm((Coor(2,:)-Coor(1,:))); 
KeL = E*A*[1,-1;-1,1]/L;
end

function [] = Plotter(Nodes,Elements,Types,Sections,Materials,dofs,P)
 K = asmkm(Nodes,Elements,Types,Sections,Materials,dofs);
 U  = K\P;
ForcesLCS = elemforces(Nodes,Elements,Types,Sections,Materials,dofs,U,[]);
f1 = figure;
set(f1, 'DefaultFigureRenderer', 'painters');
axis equal; colormap winter; axis off;
hold on
for i = 1:size(Elements,1)
    x1 = Nodes(Elements(i,5),2); x2 = Nodes(Elements(i,6),2);
    y1 = Nodes(Elements(i,5),3); y2 = Nodes(Elements(i,6),3);
    z1 = Nodes(Elements(i,5),4); z2 = Nodes(Elements(i,6),4);
   
    
FN = ForcesLCS(i,1);
Fmin = -5.6888e+04;  Fmax = 5.6888e+04;
   RR = (FN<=0)+(FN>0 )*(1-FN/Fmax);
   GG = (FN<0 )*(1-FN/Fmin)+(FN==0)+(FN>0 )*(1-FN/Fmax);
   BB = (FN>=0)+(FN<0 )*(1-FN/Fmin);
   line([x1,x2],[y1,y2],[z1,z2],'Color',[RR GG BB],'LineWidth',1);
end
view([38 30]);
hold off;
set(f1, 'DefaultFigureRenderer', 'painters');
end
