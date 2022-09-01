function masterFile(chapter)
% Reproduce results from chapters 3 to 7 with the command:
% masterFile("3")
% masterFile("4")
% masterFile("5")
% masterFile("6a")
% masterFile("6b")
% masterFile("7")

addpath(genpath('./basis functions/'));
addpath(genpath('./stabil/'));
addpath(genpath('./stabilWD/'));







if strcmp(chapter,"3")
    addpath(strcat('./chapter',{' '},chapter,'/'))
    testCaseMorley(4);
    CantileverTest();
    testCantileverPlate();
end

if strcmp(chapter,"4")
    addpath(strcat('./chapter',{' '},chapter,'/'));
    optimalDomeFmin(8);
end

if ~isempty(intersect(chapter,["5","6a","6b","7"]))
    addpath(strcat('./chapter',{' '},chapter,'/'))
    Mfunc = @ModelFunc;

    parameters                  = struct;
    parameters.derivatives      = true; 
    parameters.algorithm        = 'fmincon';
    parameters.Obj              = {'volume'}; 
    parameters.plot             = true;
    parameters.regularization   = false;
    parameters.topologyChange   = false;
    parameters.merge            = false;
    parameters.cMax             = false;
end

if strcmp(chapter,"5")
    parameters.maxIter          = 1000;
    parameters.Cons             = {'compliance'};
    parameters.topologyChange   = true;
    parameters.cMax             = true;
    parameters.cMaxV            = 30;
    penalMat = [1;1;1;1;1;2;3];
    jointCostMat = [-1;-0.5;0;0.5;1;0;0];
	for scenario = 1:size(penalMat(:),1)
        results = struct; results.parameters = parameters;
        nProb = 6;
        for problemNr = 1:nProb
            warning('off','all')
            load(strcat('./chapter',{' '},chapter,'/infochapter5'),'info')
            load(strcat('./chapter',{' '},chapter,'/Testcase',num2str(problemNr)),'Model');
            warning('on','all')
            Lfunc = @(N,E,T,P,varargin) LoadVoronoi(N,E,P,varargin{:});
            Model   = Mfunc('changeLoad',Model,info);
            Model.info.power     = penalMat(scenario);
            Model.info.jointCost = jointCostMat(scenario);
            Model.info.shift     = false;
            Problem = Setup(Model);
            results.(strcat('case',num2str(problemNr))) = optimization(Model,Problem,Mfunc,Lfunc,parameters);
        end
        Plotter(results,Mfunc,Lfunc,nProb);
	end
end

if strcmp(chapter,"6a")
    parameters.maxIter          = 200;
    parameters.Cons             = {'displacements','stress','buckling'};
    parameters.regularization   = true;
    stepMat = [1e-5;1.5];
    results = struct; results.parameters = parameters;
	for scenario = 1:size(stepMat(:),1)
		warning('off','all')
        load(strcat('./chapter',{' '},chapter,'/Testcase',num2str(1)),'Model');
        warning('on','all')
        Lfunc                = @(N,E,T,P,varargin) LoadVoronoi(N,E,P,varargin{:});
        Model                = ModelFunc('init',Model);
        Model.info           = struct;
        Model.info.power     = 1;
        Model.info.jointCost = 0;
        Model.info.nClusterR = 10;
        Model.info.nClusterS = 40;
        Model.info.Rad       = 34; 
        Model.info.shift     = true;
        Model.info.stepP     = stepMat(scenario);
        [Problem] = Setup(Model);
        results.(strcat('case',num2str(scenario))) = optimization(Model,Problem,Mfunc,Lfunc,parameters);
	end
    Fex = [-681, 681; -17, 17];
    Plotter(results,Mfunc,Lfunc,size(stepMat(:),1),Fex);
end

if strcmp(chapter,"6b")
    parameters.maxIter          = 200;
    parameters.Cons             = {'displacements','stress','buckling','Hmax'};
    stepMat = [1e-3;5];
    results = struct; results.parameters = parameters;
	for scenario = 1:size(stepMat(:),1)
        warning('off','all')
        load(strcat('./chapter',{' '},chapter,'/Testcase',num2str(1)),'Model');
        warning('on','all')
        Lfunc                = @(N,E,T,P,varargin) LoadVoronoi(N,E,P,varargin{:});
        Model.info.power     = 1;
        Model.info.jointCost = 0;
        Model.info.nClusterR = 10;
        Model.info.nClusterS = 40;
        Model.info.shift     = false;
        Model.info.stepP     = stepMat(scenario); 
        Problem = Setup(Model);
        results.(strcat('case',num2str(scenario))) = optimization(Model,Problem,Mfunc,Lfunc,parameters);
	end
    Fex = [-681, 681; -17, 17];
    Plotter(results,Mfunc,Lfunc,size(stepMat(:),1),Fex);
end

if strcmp(chapter,"7")
    parameters.maxIter          = 200;
    parameters.Cons             = {'displacements','stress','buckling','surfaceFunc','deflectionCladding','compliance'};
    parameters.regularization   = true;
    parameters.merge            = true;
    parameters.Rmerge           = 0.01;
    deltaMat = [0.1;0.1;0.1;0.5;0.5;0.5];
    LC2Mat   = [0;0.5;1;0;0.5;1];
	for scenario = 1:size(deltaMat(:),1)
        results = struct; results.parameters = parameters;
        nProb = 6;
        for problemNr = 1:nProb
            warning('off','all')
			load(strcat('./chapter',{' '},chapter,'/Testcase',num2str(problemNr)),'Model');
            warning('on','all')
            Lfunc                = @(N,E,T,P,varargin) LoadShape(N,T,P,varargin{:});
            Model.info           = struct;
            Model.info.power     = 1;
            Model.info.jointCost = 0;
            Model.info.nClusterR = 10;
            Model.info.nClusterS = 10;
            Model.info.delta     = deltaMat(scenario);
            Model.info.LC2       = LC2Mat(scenario);  
            Model.info.shift     = false;
            Model.info.Rad       = 34;  
            [Problem] = Setup(Model);
            Model = Mfunc('postMerge',Model.Node,Model.Elem,Model.Trig,Problem.x,Model.info);
            [Problem] = Setup(Model);
            results.(strcat('case',num2str(problemNr))) = optimization(Model,Problem,Mfunc,Lfunc,parameters);
        end
        Fex = [-920, 920; -180, 180];
        Plotter(results,Mfunc,Lfunc,nProb,Fex);   
	end
end


if strcmp(chapter,"Experiment")
    pathname = strcat('./chapter',{' '},num2str(5),'/');
    addpath(pathname{1})
    Mfunc = @ModelFunc;
    parameters                  = struct;
    parameters.derivatives      = true; 
    parameters.algorithm        = 'fmincon';
    parameters.Obj              = {'volume'}; 
    parameters.plot             = true;
    parameters.regularization   = false;
    parameters.topologyChange   = false;
    parameters.merge            = false;
    parameters.cMax             = false;
    parameters.maxIter          = 1000;
    parameters.Cons             = {'compliance'};
    nVerts = 6;
    alp = pi/nVerts;
    info.domain = polyshape([0,0;cos(alp),0;cos(alp),sin(alp)]*28);
    info.lMaxF    = 15;
    info.nStep    = 4;
    verts         = info.domain.Vertices;
    info.refNode  = [0,0,-1e10];
    info.suppNode = {[verts(2,:),0],3};
    info.suppRule = {@(nod) nod(abs(nod(:,4))<1e-5,1), [0.01,0.02,0.03]};


    info.loadNode = {[verts(2,:),0],[0,0,0,0,0,0],1};
    info.loadRule = {   @(nod) nod(1:end-1,1),[0,0,-3.92,0,0,0],1;
                        @(nod) nod(1:end-1,1),[0,-1.68,0,0,0,0],2}; 
    results = struct; results.parameters = parameters;
    Lfunc = @(N,E,T,P,varargin) LoadVoronoi(N,E,P,varargin{:});
    unsuccessful = true;
    while unsuccessful
        try
            Model   = Mfunc('init',info);
            unsuccessful = false;
        catch
        end
    end
    Model.info.power     = 1;
    Model.info.jointCost = 0;
    Problem = Setup(Model);
    results.(strcat('case',num2str(1))) = optimization(Model,Problem,Mfunc,Lfunc,parameters);
    Plotter(results,Mfunc,Lfunc,nProb);
end

