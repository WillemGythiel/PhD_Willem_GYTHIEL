function Plotter(results,Mfunc,Lfunc,nProblems,varargin)
perc = cell(1,nProblems);
fValAll = zeros(1010,nProblems);
for problemNr = 1:nProblems
    res = results.(strcat('case',num2str(problemNr)));
    Model = res.ModelFin; Problem = res.ProblemFin; fVal = res.fval; x = res.ProblemFin.x;
    if nargin > 4
       Fex =  varargin{1};
       plotDevice1(Model,Problem,x,Fex,Mfunc,Lfunc)
    end
    perc{problemNr} = plotDevice2(Model,Problem,x,Mfunc,Lfunc);
    Weight = sum(elemvolumes(Model.Nodes,Model.Elements,Model.Types,Model.Sections),1)*7.85;
%     HarryPlotter2(Model.Nodes,Model.Elements,Model.Sections)
    fValAll(1:size(fVal,1),problemNr) = fVal*Weight/fVal(end);                   
end
fValAll = fValAll(sum(fValAll,2)~=0,:);

figure
hold on
for problemNr = 1:nProblems
    plot((1:size(fValAll,1))',fValAll(:,problemNr))
end
hold off
xlabel('iteration') 
ylabel('Weight') 

figure
perc = cell2mat(perc)';
bar(perc(1:size(Model.Pdis,2):end,:),'stacked')
xlabel('case') 
ylabel('compliance contribution') 
end



function [] = plotDevice1(Model,Problem,x,Fex,Mfunc,Lfunc)
M = Mfunc('update',Model,Problem,x);
K = asmkm(M.Nodes,M.Elements,M.Types,M.Sections,M.Materials,M.dofs);
P = Lfunc(M.Nodes,M.Elements,M.Trigs,M.Pdis);
nlc = size(M.Pdis,2);       nDof = size(M.dofs,1); 
DOFs = round(6*(round(M.dofs)-1)+100*mod(M.dofs,1));
P = reshape(P(DOFs,:,:),nDof,nlc);
U  =K\P;
ForcesLCS = elemforces(M.Nodes,M.Elements,M.Types,M.Sections,M.Materials,M.dofs,U);

for LC = 1:nlc
    f1 = figure();
    SubPlot1_theForce(M,ForcesLCS(:,:,LC),f1,Fex)
    f2 = figure();
    SubPlot2_justMomentsAway(M,ForcesLCS(:,:,LC),f2,Fex)
end
end

function [] = SubPlot1_theForce(M,ForcesLCS,f1,Fex)
set(f1, 'DefaultFigureRenderer', 'painters');
axis equal; colormap winter; axis off;
hold on
for i = 1:size(M.Elements,1)
    x1 = M.Nodes(M.Elements(i,5),2); x2 = M.Nodes(M.Elements(i,6),2);
    y1 = M.Nodes(M.Elements(i,5),3); y2 = M.Nodes(M.Elements(i,6),3);
    z1 = M.Nodes(M.Elements(i,5),4); z2 = M.Nodes(M.Elements(i,6),4);
FN = ForcesLCS(i,1);
Fmax = max(max(abs(FN(:))),Fex(1,2));
Fmin = -Fmax;
   RR = (FN<=0)+(FN>0 )*(1-FN/Fmax);
   GG = (FN<0 )*(1-FN/Fmin)+(FN==0)+(FN>0 )*(1-FN/Fmax);
   BB = (FN>=0)+(FN<0 )*(1-FN/Fmin);
   line([x1,x2],[y1,y2],[z1,z2],'Color',[RR GG BB],'LineWidth',1);
end
view(2)
hold off;
set(f1, 'DefaultFigureRenderer', 'painters');
end

function [] = SubPlot2_justMomentsAway(M,ForcesLCS,f1,Fex)
nS = 10;
inp = @(c1,c2,v) [c1*(1-v/nS)+(v/nS)*c2, c1*(1-(v+1)/nS)+(v+1)/nS*c2];
set(f1, 'DefaultFigureRenderer', 'painters');
axis equal; colormap winter; axis off;
Mom = ForcesLCS(:,[6,12]);
hold on
for i = 1:size(M.Elements,1)
    x1 = M.Nodes(M.Elements(i,5),2); x2 = M.Nodes(M.Elements(i,6),2);
    y1 = M.Nodes(M.Elements(i,5),3); y2 = M.Nodes(M.Elements(i,6),3);
    z1 = M.Nodes(M.Elements(i,5),4); z2 = M.Nodes(M.Elements(i,6),4);
Fmax = max(max(abs(Mom(:))),Fex(2,2));
Fmin = -Fmax;
   RR = @(Fi)(Fi<=0)+(Fi>0 )*(1-Fi/Fmax);
   GG = @(Fi)(Fi<0 )*(1-Fi/Fmin)+(Fi==0)+(Fi>0 )*(1-Fi/Fmax);
   BB = @(Fi)(Fi>=0)+(Fi<0 )*(1-Fi/Fmin);
   for j = 0:nS-1
       R1 = RR(Mom(i,1));   R2= RR(Mom(i,2));
       G1 = GG(Mom(i,1));   G2= GG(Mom(i,2));
       B1 = BB(Mom(i,1));   B2= BB(Mom(i,2));
       RGB = [mean(inp(R1,R2,j)),mean(inp(G1,G2,j)),mean(inp(B1,B2,j))];
       line(inp(x1,x2,j),inp(y1,y2,j),inp(z1,z2,j),'Color',RGB,'LineWidth',1);
   end
end
view(2)
hold off;
set(f1, 'DefaultFigureRenderer', 'painters');
end

function [res] = plotDevice2(Model,Problem,x,Mfunc,Lfunc)
M = Mfunc('update',Model,Problem,x);
K = asmkm(M.Nodes,M.Elements,M.Types,M.Sections,M.Materials,M.dofs);
Kinv = K\speye(size(K));
P = Lfunc(M.Nodes,M.Elements,M.Trigs,M.Pdis);
nlc = size(M.Pdis,2);       nDof = size(M.dofs,1); 
DOFs = round(6*(round(M.dofs)-1)+100*mod(M.dofs,1));
P   = reshape(P(DOFs,:,:),nDof,nlc);
Uall = Kinv*P;
M.Sections(:,3:4) = 1e10;
SectionsA  = M.Sections*diag([1,1,0,0,0,0,0,0,0,0,0]);
SectionsIx = M.Sections*diag([1,1e-10,0,0,1,0,0,0,0,0,0]); 
SectionsIy = M.Sections*diag([1,1e-10,1,1,0,1,0,1,1,1,1]); SectionsIy(:,3:4) = Inf;
SectionsIz = M.Sections*diag([1,1e-10,1,1,0,0,1,1,1,1,1]); SectionsIz(:,3:4) = Inf;
KA = asmkm(M.Nodes,M.Elements,M.Types,SectionsA,M.Materials,M.dofs);
KIx = asmkm(M.Nodes,M.Elements,M.Types,SectionsIx,M.Materials,M.dofs);
KIy = asmkm(M.Nodes,M.Elements,M.Types,SectionsIy,M.Materials,M.dofs);
KIz = asmkm(M.Nodes,M.Elements,M.Types,SectionsIz,M.Materials,M.dofs);

res = zeros(4,size(Uall,2));

for lc = 1:size(Uall,2)
    U = Uall(:,lc);
    comp  = U'*K*U;
    conA  = U'*KA*U/comp;
    conIx = U'*KIx*U/comp;
    conIy = U'*KIy*U/comp;
    conIz = U'*KIz*U/comp;
    res(:,lc) = [conA;conIx;conIy;conIz];
end
end

% function [Lav] = Lcal(M)
% nelm = size(M.Elements,1);
% L = 0;
% for i = 1:nelm
%     x1 = M.Nodes(M.Elements(i,5),2); x2 = M.Nodes(M.Elements(i,6),2);
%     y1 = M.Nodes(M.Elements(i,5),3); y2 = M.Nodes(M.Elements(i,6),3);
%     z1 = M.Nodes(M.Elements(i,5),4); z2 = M.Nodes(M.Elements(i,6),4);
%    
%     L = L+((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)^0.5;
% end
% Lav = L/nelm;
% end
