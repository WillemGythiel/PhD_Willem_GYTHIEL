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
optimValues.funccount = 0;
optimValues.fval = 0;
state = 'init';
% history = struct;
while true
  state = 'interrupt';
  [f0,df0]        = fun(x);
  [f,feq,df,dfeq] = nonlcon(x);
  if ~(isempty(feq)&&isempty(dfeq)&&isempty(A)&&isempty(b)&&...
                                    isempty(Aeq)&&isempty(beq))
      disp('Invalid input');
      break
  end
  xmin = max(lb,x-options.StepSize);
  xmax = min(ub,x+options.StepSize);
  [xnew,~,~,~,mp,~] = mma(x,xmin,xmax,f0,f,df0',df',mp,a0,a,c,d,[],[],[],s0,si,sd);
  change = max(abs(xnew-x));
  optimValues.fval = f0;
  if options.display == 'iter'
    disp([' It.: ' sprintf('%4i',optimValues.iteration) ...
          ' ch.: ' sprintf('%6.3f',change )...
          '      ' sprintf(datestr(now,'HH:MM:SS'))]);
  end
  if  ~isempty(options.OutputFcn)
    options.OutputFcn(x,optimValues,state);
  end

  state = 'iter';
  if change>=options.StepTolerance && optimValues.iteration+1<options.MaxIterations
    optimValues.iteration = optimValues.iteration+1;
    optimValues.funccount = optimValues.funccount+1;
    x = xnew;
  else
    varargout{1} = f0;
    varargout{2} = change<options.StepTolerance*4;
    output = optimValues;
    output.iterations = optimValues.iteration;
    varargout{3} = output;
    break;
  end
end
 state = 'done';
 
end