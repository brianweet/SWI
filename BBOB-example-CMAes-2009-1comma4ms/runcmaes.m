function [XMIN, FMIN, COUNTEVAL, STOPFLAG, OUT, BESTEVER] = ...
            runcmaes(FUN, DIM, ftarget, maxfunevals)

  clear opts;
  
  % set unmirrored CMA without serialism although standard settings:
  opts.Serial = 'on';
	opts.Mirrored = 'on';
	opts.PopSize=4;
  opts.ParentNumber=1;
  opts.evalpar = 0;
  

  opts.restarts = 0; % no restarts with increasing popsize
  %opts.maxiter = maxfunevals; % allow full number of function evaluations
	opts.maxiter = 2 * (100 + ceil(1000*DIM*sqrt(DIM)));

  % ultimate termination option
  opts.maxfunevals = maxfunevals;

  % speed and output options
  opts.stopfitness = ftarget;

  opts.savevar = 'off';
  opts.dispfinal = 0;
  opts.dispmod = 0;
  opts.logmod = 0;
  
  % use parameters from Niko's BiPOP-CMAES from BBOB'2009:
	opts.TolHistFun = '1e-12';
  
  % restarts as long as we have function evaluations left:
  FMIN=inf;  
  while(FMIN > ftarget) && (opts.maxfunevals > 1)
  
      [XMIN, FMIN, COUNTEVAL, STOPFLAG, OUT, BESTEVER] = ...
        cmaes(FUN, ['8*rand(' num2str(DIM) ', 1) - 4'], 2, opts);
  
    opts.maxfunevals=opts.maxfunevals - COUNTEVAL;
  end

