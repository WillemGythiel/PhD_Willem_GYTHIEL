function [x,varargout] = MMA_func(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,varargin)
x    = x0; 
if nargin<10
    options = optimoptions_MMA('MMA');
else
    options = varargin{1};
end
% MMA PARAMETERS
a0      = 1;                            % Weight for artificial variable z in obj function
a       = 0;                            % Weight for artificial variable z in constraints
c       = 10000;                        % Constraint relaxing parameter
d       = 0;                            % Parameter for least squares problems
mp      = [];                           % Internal MMA parameters
s0      = 0.01;                         % Initial asymptote location
si      = 1.01;                         % Asymptote increase factor
sd      = 0.7;                          % Asymptote decrease factor

% OPTIMIZATION LOOP
optimValues = struct;
optimValues.iteration = 0;
optimValues.funcCount = 0;

if strcmp(options.display,'iter')
    disp("   Iter   fCount    F(x)         conViol       step    ")
end
% state = 'init';
% history = struct;
while true
  % state = 'interrupt';
  state = 'iter';
  [f0,df0]        = fun(x);
  [f,feq,df,dfeq] = nonlcon(x);
  if ~(isempty(feq)&&isempty(dfeq)&&isempty(A)&&isempty(b)&&...
                                    isempty(Aeq)&&isempty(beq))
      disp('Invalid input');
      break
  end
  optimValues.constrviolation = max(max(f(:)),0);
  xmin = max(lb,x-options.StepSize);
  xmax = min(ub,x+options.StepSize);
  
  [xnew,~,~,~,mp,~] = mma(x,xmin,xmax,f0,f,df0',df',mp,a0,a,c,d,[],[],[],s0,si,sd);
  change = max(abs(xnew-x));
  optimValues.fval = f0;
  if strcmp(options.display,'iter')
      printM = [optimValues.iteration,optimValues.funcCount,f0,optimValues.constrviolation,change];
      fprintf('%6i %6i %12.6f %12.3f %12.4f\n', printM);
  end
  
  stop = false;
  if  ~isempty(options.OutputFcn)
    stop = options.OutputFcn(x,optimValues,state);
  end

  optimValues.iteration = optimValues.iteration+1;
  optimValues.funcCount = optimValues.funcCount+1;
  x = xnew;  
  
  if stop || optimValues.iteration>options.MaxIterations || change<=options.StepTolerance
  output = struct;
  output.iteration = optimValues.iteration;
  output.funcCount = optimValues.funcCount;
  output.algorithm = 'MMA';
  output.constrviolation = max(max(f(:)),0);
  output.stepsize = options.StepSize;
      if stop
        exitflag = -1;
      elseif optimValues.iteration>options.MaxIterations
        exitflag = 0;
      elseif change<=options.StepTolerance
        exitflag = 1;
      end
  varargout{1} = f0;
  varargout{2} = exitflag;
  varargout{3} = output;
  break;
  end
end
% state = 'done';
 
end