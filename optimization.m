function [history] = optimization(Model,Problem,Mfunc,Lfunc,parameters)
if parameters.plot
    h = figure; set(gcf,'color','w');   axis tight manual  
end

history = struct; history.x = []; history.fval = []; history.cons = []; history.iter = 0; history.time = 0;


Vmax = objCon(Problem.x,'objective',{'volume'},Mfunc,Lfunc,Model,Problem,0,0);
cCmaxC = 5*max(objCon(Problem.x,'objective',{'complianceCladding'},Mfunc,Lfunc,Model,Problem,0,0));
if parameters.cMax
    Cmax = parameters.cMaxV;
else
    Cmax = 1.1*objCon(Problem.x,'objective',{'compliance'},Mfunc,Lfunc,Model,Problem,0,0);
end

flag = -1; % -1 means the process was interrupted by the outpunt function for merging, which is fine
while flag == -1
    options = makeOptions(parameters);
    Problem.Vmax = Vmax;    Problem.cCmaxC = cCmaxC;    Problem.Cmax = Cmax;
    fun     = @(x) objCon(x,'objective', parameters.Obj,Mfunc,Lfunc,Model,Problem,parameters.derivatives);
    nonlcon = @(x) objCon(x,'constraint',parameters.Cons,Mfunc,Lfunc,Model,Problem,parameters.derivatives,parameters.regularization);
    tic
    if strcmp(parameters.algorithm,'fmincon')
        [xNew,~,flag,~] = fmincon(fun,Problem.x,[],[],[],[],Problem.lb,Problem.ub,nonlcon,options);  
    elseif strcmp(parameters.algorithm,'ga')
        [xNew,~,flag,~] = ga(fun,size(Problem.x,1),[],[],[],[],Problem.lb,Problem.ub,nonlcon,options);
    elseif strcmp(parameters.algorithm,'MMA')
        [xNew,~,flag,~] = MMA_func(fun,Problem.x,[],[],[],[],Problem.lb,Problem.ub,nonlcon,options);
    end
    Problem.x = xNew(:);
    history.time = history.time+toc;
    
    

    if parameters.topologyChange
        nelmOld = size(Model.Elements,1);
        [Model,Problem] = topology_changer(Mfunc,Model,Problem);
        nelmNew = size(Model.Elements,1);
        if nelmNew<nelmOld && history.iter < parameters.maxIter
            flag = -1;
        end
    end    
    
    if parameters.merge
      [Model,Problem] = MergeProcess(Mfunc,Model,Problem,parameters.Rmerge);
    end
end

function [state,options,optchanged] = updateStructGA(options,state,flag)
optchanged = false;
switch flag
    case 'iter'
       % Concatenate current point and objective function
       % value with history. x must be a row vector.
       aux = (1:numel(state.Score));
       ind = min(aux(state.Score==state.Best(end)));
       history.fval = [history.fval; state.Best(end)];
       history.x    = [history.x; state.Population(ind,:)];
       history.iter = history.iter+1;
       if parameters.plot
           filename = 'Gridshell.gif';
           M = GSModel('update',Model,Problem,x);
           MakeGIF_GS(M.Nodes,M.Elements,M.Trigs,M.Sections,h,filename,history.iter)
       end
   otherwise
end
end

function [stop] = updateStruct(x,optimValues,state)
stop = false;
switch state
    case 'iter'
       % Concatenate current point and objective function
       % value with history. x must be a row vector.
       history.fval = [history.fval; optimValues.fval];
       history.x    = [history.x; x];
       history.cons = [history.cons; optimValues.constrviolation];
       history.iter = history.iter+1;
       if parameters.plot
           filename = 'Gridshell.gif';
           M = Mfunc('update',Model,Problem,x);
           MakeGIF_GS(M.Nodes,M.Elements,M.Trigs,M.Sections,h,filename,history.iter)
       end
       if parameters.merge
       	stop = merge_checker(Mfunc,Model,Problem,parameters.Rmerge,x);
       end
   otherwise
end
end

function options = makeOptions(parameters)
if strcmp(parameters.algorithm,'fmincon')
    options = optimoptions('fmincon','OutputFcn',@(x,o,s) updateStruct(x,o,s),...
    'SpecifyObjectiveGradient',parameters.derivatives,'SpecifyConstraintGradient',parameters.derivatives,...
    'FunctionTolerance',1e-8,'StepTolerance',1e-12,...
    'MaxIterations',parameters.maxIter-history.iter,'display','iter'); %,'PlotFcn', {@optimplotfval}
elseif strcmp(parameters.algorithm,'ga')
    options = optimoptions('ga','OutputFcn',@(o,s,f) updateStructGA(o,s,f),...
    'FunctionTolerance',1e-6,...
    'display','iter');
elseif strcmp(parameters.algorithm,'MMA')
    options = optimoptions_MMA('MMA','OutputFcn',@(x,o,s) updateStruct(x,o,s),...
    'FunctionTolerance',1e-6,'StepTolerance',1e-6,...
    'MaxIterations',parameters.maxIter-history.iter,'display','iter');
end
end


history.ModelFin = Model;
history.ProblemFin = Problem;
end




function MakeGIF_GS(Nodes,Elements,Trigs,Sections,h,filename,iter)
PlotStructure(Nodes,Elements,Trigs,Sections,['Iteration: ' num2str(iter)])
% Capture the plot as an image 
frame = getframe(h); 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256); 
% Write to the GIF File 
n = iter;
if n == 1 
  imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0); 
else 
  imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0); 
end 
end

function PlotStructure(Nodes,Elements,Trigs,Sections,varargin)
Bars = Elements(:,[1,5,6]);
xSec = 60*Sections(Elements(:,3),2);
% figure
trisurf(Trigs(:,2:4), Nodes(:,2), Nodes(:,3), Nodes(:,4), 'visible', 'off','facecolor',[0.5,0.5,0.5]);
axis equal; colormap winter; axis off; 
% view(2);
view(45,45)
if nargin > 4
    TitleText = varargin{1};
    title(TitleText)
end
ax = gca;
ax.TitleHorizontalAlignment = 'left';
hold on
for i = 1:size(Bars,1)
    x1 = Nodes(Bars(i,2),2); x2 = Nodes(Bars(i,3),2);
    y1 = Nodes(Bars(i,2),3); y2 = Nodes(Bars(i,3),3);
    z1 = Nodes(Bars(i,2),4); z2 = Nodes(Bars(i,3),4);
   line([x1,x2],[y1,y2],[z1,z2],'LineWidth',xSec(i),'color',[0,0,0]); 
end
hold off
drawnow
end

function [Model,Problem] = MergeProcess(Mfunc,Model,Problem,varargin)
if nargin <4
    Dmin = 0.1;
else
    Dmin = varargin{1};
end
if nargin <5
    [Model,Problem] = Merging(Mfunc,Model,Problem,Dmin);
else
    test = varargin{2}; tol = 0.05;
    unsuccessful = true;
    while unsuccessful
        try
            [ModelM,ProblemM] = Merging(Mfunc,Model,Problem,Dmin);
            assert(abs((test(Problem.x)-test(ProblemM.x))/test(Problem.x))<tol);
            unsuccessful = false;
        catch
            Dmin = Dmin/2;
        end
    end
    Model = ModelM; Problem = ProblemM;
end
end

function [Model,Problem] = Merging(Mfunc,Model,Problem,Dmin)
Model = Mfunc('update',Model,Problem,Problem.x);
[Vert,Peri] = Model.cache.get_protected_nodes(Model.Node);
[Node,Elem,Trig,x] = mergerGS(Model.Node,Model.Elem,Model.Trig,Problem.x(end-size(Model.Sections,1)+1:end),Dmin,Vert,Peri);
Model = Mfunc('postMerge',Node,Elem,Trig,x,Model.info);
Problem = Setup(Model);
end


function [stop4merge] = merge_checker(Mfunc,Model,Problem,Dmin,x)
Model = Mfunc('update',Model,Problem,x);
[x,y] = ndgrid(Model.Node(:,2), Model.Node(:,3));
dij = ((x-x').^2+(y-y').^2).^0.5;
stop4merge = sum(nonzeros(dij)<Dmin,1)>0;
end


function [Model,Problem] = topology_changer(Mfunc,Model,Problem)
Model = Mfunc('update',Model,Problem,Problem.x);
Model = Mfunc('topology',Model,Problem,Problem.x);
Problem = Setup(Model);
end



