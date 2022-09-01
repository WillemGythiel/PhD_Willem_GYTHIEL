function [] = testCaseMorley(i)
for w = 1:i
%% INPUT PROPERTIES
n = 2^(w+1);
La = 100; Lb = La;    theta = pi/6;
t  = 0.1;   E  = 1.09e6;    nu = 0.3;   Load = 1;
D  = E*t^3/(12*(1-nu^2));
%% DEFINITION OF THE NODES
num  = (n+1)^2;    
    nID = (1:num)';                                                        % nodal IDs
    nx  = kron(ones(n+1,1),(0:La/n:La)')+kron((0:n)',Lb/n*cos(theta)*ones(n+1,1));                                  % nodal coordinates
    ny  = kron((0:Lb/n:Lb)',ones(n+1,1))*sin(theta);
    nz  = zeros(num,1);
Nodes = [nID,nx,ny,nz];                                                    % difinition of the nodes
%% DEFINITION OF THE ELEMENTS
inds=(1:n^2+n)'; trigV = inds(mod(inds,n+1)~=0);                           % triagonal indices
    upTrigs   = horzcat(trigV,   trigV+1,    trigV+n+1);
    downTrigs = horzcat(trigV+1, trigV+n+2,  trigV+n+1);
Trigs = vertcat(upTrigs, downTrigs);

%% DEFINITION OF THE DOFs
rv1  = ny == 0;                      rv3  = mod(nID,n+1) == 1;  
rv2  = ny == Lb*sin(theta);          rv4  = mod(nID,n+1) == 0;
Edge = (rv1 | rv2 | rv3 | rv4);      Crnr = (rv1 | rv2)&(rv3 | rv4);
%     dofz = 3*(nID(Edge==0)-1)+1;
    dofz = 3*(nID(Crnr==0)-1)+1;
    dofr = Addkron(3*(nID-1),(2:3)');
dofs = vertcat(dofz, dofr);
%% DEFINITION OF THE LOADS
vw1 = ny == 0.5*Lb*sin(theta);      vw2 = mod(nID,n+1) == n/2+1;
vwP = kron([1;0;0],vw1 & vw2);
Pload = Load*(1-Edge/2-Crnr/4)*La*Lb*sin(theta)/(n^2);
nDOF = 3*size(Nodes,1);      nTRI = size(Trigs,1);
P = sparse(3*(nID-1)+1,0*nID+1,Pload,nDOF,1);
%% DEFINITION OF THE STIFFNESS MATRIX
  Mi = zeros(81*nTRI,1);      Mj = Mi;  Mv = Mi;                                  % definition of stiffness matrix K
for q = 1:size(Trigs,1)
    C = (Nodes(Trigs(q,:),2:4));
    C = C-repmat(C(1,:),3,1);
    [Ke] = KeLCS_TrigPlateDKT(C, t, E, nu);
    pos  = Addkron(3*(Trigs(q,:)-1),(1:3));
    ind  = 81*(q-1)+(1:81);
    Mi(ind) = kron(ones(1,9),pos);
    Mj(ind) = kron(pos,ones(1,9));
    Mv(ind) = Ke(:);
end
K = sparse(Mi,Mj,Mv,nDOF,nDOF);
U = sparse(nDOF,1);
U(dofs) = K(dofs,dofs)\P(dofs);
wMax = U(3*(nID(vwP==1)-1)+1);
wTest = Load*La^4/D;
wFin  = wMax/wTest;
% figure
%     trisurf(Trigs,nx,ny,nz+0.2*La/wMax*U(1:3:end));
%     axis equal; colormap winter; axis off;
  history.subdiv(:,w) = n;
  history.wFin(:,w) = wFin;
end
graphMaker(history.subdiv,history.wFin,'Morley skew plate - DKT',...
           'subdivision',{'Dimensionless';'mid deflection'},2);
end

function [] =graphMaker(xval,yval,titleT,labelX,labelY,opt)
% This function makes simple, clear line graphs, inspired bij J-L Doumont. 
%
% Input:    xval   = values on x axis       (vector of doubles)
%           yval   = correspondig y values  (vector of doubles)
%           titleT = the title              (string or cell of strings)
%           labelX = label for x axis       (string or cell of strings)
%           labelY = label for y axis       (string or cell of strings)
%           opt    = 1: y axis shows range of values
%                    2: y axis shows lowest, highest and trend value
%
%--------------------------------SIDENOTE---------------------------------%
% The y and x axies are drawn separately in 2 different subplots, which
% allows displaying a gap between the axes if they do not intersect at
% the origin (i.e. {0,0}). The data is plotted next to the y axis in the
% upper plot. The axes are aligned and scaled to nicely fill the figure.
% The overall layout of the grapgh is minimalistic and clear (white
% background, dicrete gray axes...).
%-------------------------------------------------------------------------%
figure
%% input
linewidth  = 1.2;
fontsize   = 12;
fontname   = 'arial';
fontweight = 'normal';
axCol      = [0 0 0];
rounding   = 10;
%% upperplot
% plot
ax1 = subplot(2,1,1);
% x increment
 incX = (max(xval)-min(xval))/5;
 tenX = floor(log10(incX));
 fivX = max([round(fix(incX/(10^tenX))/5)*5,...
             round(fix(incX/(10^tenX))/2)*2*...
            (round(fix(incX/(10^tenX))/2)*2<4),1]);
xstep = fivX*10^tenX;
% y increment
 incY = (max(yval)-min(yval))/5;
 tenY = floor(log10(incY));
 fivY = max([round(fix(incY/(10^tenY))/5)*5,...
             round(fix(incY/(10^tenY))/2)*2*...
            (round(fix(incY/(10^tenY))/2)*2<4),1]);
ystep = fivY*10^tenY;
% plots
hold on
yTrend = 0.000408;
line([min(xval) max(xval)],[yTrend    yTrend   ],...
    'Color',[1,0.75,0],'LineWidth',0.8*linewidth);
plot(xval,yval,'LineWidth',2*linewidth);
hold off
% adjusting the settings
ax = gca; ax.Box = 'off'; ax.Layer = 'top'; set(gcf,'color','w');
ax.XAxis.Color = 'None'; 
ax.YAxis.Color = axCol;
ax.LineWidth   = linewidth;
% font
ax.FontSize   = fontsize;
ax.FontName   = fontname;
ax.FontWeight = fontweight;
% markers
ax.TickDir = 'out';

if opt ==1
% limits
 yAxMin = 0+floor(min(yval)/ystep)*ystep;
 yAxMax = 0+ ceil(max(yval)/ystep)*ystep;
% markers
ax.YTick = yAxMin:ystep:yAxMax;
  xAxMin = 0+floor(min(xval)/xstep)*xstep;
  xAxMax = 0+ ceil(max(xval)/xstep)*xstep;
 ax.XLim = [xAxMin xAxMax];
 ax.YLim = [yAxMin yAxMax];
ax.XTick = xAxMin:xstep:xAxMax;
else

yAxMin = 0+round(min(min(yval),yTrend),rounding);    
yAxMax = 0+round(max(max(yval),yTrend),rounding);
% markers
ax.YTick = unique(round([yAxMin,yval,yTrend,yAxMax],rounding)); 
  xAxMin = 0+min(xval);
  xAxMax = 0+max(xval);
 ax.XLim = [xAxMin xAxMax];
 ax.YLim = [yAxMin yAxMax];
ax.XTick = unique(xval);
end


% title
t = title(titleT);
% label
y = ylabel(labelY,'FontSize',fontsize,'FontWeight',fontweight,'Color',axCol, 'Rotation', 0);
% putting it all in the right place
tpos = get(t, 'Position');
ypos = get(y, 'Position');
ydim = get(y, 'Extent');
set(y, 'Position', [ax.XLim(1) tpos(2)+0.2*ydim(4) ypos(3)])
set(t, 'Position', [tpos(1)  tpos(2)+0.2*ydim(4) tpos(3) ])


%% lowerplot
ax2 = subplot(2,1,2);
% adjusting the settings
ax = gca; ax2.Box = 'off'; ax.Layer = 'top'; set(gcf,'color','w');
ax2.XAxis.Color = axCol; 
ax2.YAxis.Color = 'None';
ax2.LineWidth   = linewidth;
% font
ax2.FontSize  = fontsize;
ax2.FontName  = fontname;
ax2.FontWeight = fontweight;


if opt ==1
% markers
ax2.TickDir = 'out';
ax2.XTick = xAxMin:xstep:xAxMax;
% limits
ax2.XLim = [xAxMin xAxMax];
else
    ax2.TickDir = 'out';
  xAxMin = 0+min(xval);
  xAxMax = 0+max(xval);
 ax.XLim = [xAxMin xAxMax];
ax.XTick = unique(xval);
end

% labels
xlabel(labelX,'FontSize',fontsize,'FontWeight',fontweight,'Color',axCol);
%% aligning the upper and lower plot
dims = get(ax, 'Position');
difH = 1.75*dims(2);
addH = get(ax, 'OuterPosition')+get(ax2, 'OuterPosition');
totH = addH(4);
remH = totH-2*difH;
dims2 = get(ax2, 'Position');
dims(4) = remH;
dims(2) = difH+0.025*totH*(yAxMin~=0);
dims2(2) = difH-0.025*totH*(yAxMin~=0);
dims2(4) =  0;
set(ax1, 'Position', dims);
set(ax2, 'Position', dims2);

end


