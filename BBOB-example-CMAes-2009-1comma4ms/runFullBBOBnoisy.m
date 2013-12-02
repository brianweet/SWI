% runs an entire experiment for benchmarking MY_OPTIMIZER
% on the noise-free testbed. fgeneric.m and benchmarks.m
% must be in the path of Matlab/Octave
% CAPITALIZATION indicates code adaptations to be made

addpath('./');  % should point to fgeneric.m etc.

opt.algName = '(1,4_m^s)CMA-ES';
opt.comments = 'serial (1,4)CMA-ES with restarts, mirrored version';
maxfunevals = '1e4 * dim';  % SHORT EXPERIMENT, takes overall three minutes 

more off;  % in octave pagination is on by default

t0 = clock;
rand('state', sum(100 * t0));

datapath = './outputFiles1komma4mirroredserialwith1e4timesDimFunEvalNoisy';  % different folder for each experiment

for dim = [20,10,5,3,2]  % small dimensions first, for CPU reasons
  for ifun = benchmarksnoisy('FunctionIndices')
    for iinstance = 1:15 %[1:5, 1:5, 1:5]  % first 5 fct instances, three times
      fgeneric('initialize', ifun, iinstance, datapath, opt);
      
      runcmaes('fgeneric',dim,fgeneric('ftarget'),eval(maxfunevals)); 
     
      disp(sprintf(['  f%d in %d-D, instance %d: FEs=%d,' ...
                    ' fbest-ftarget=%.4e, elapsed time [h]: %.2f'], ...
                   ifun, dim, iinstance, ...
                   fgeneric('evaluations'), ...
                   fgeneric('fbest') - fgeneric('ftarget'), ...
                   etime(clock, t0)/60/60));
      fgeneric('finalize');
    end
    disp(['      date and time: ' num2str(clock, ' %.0f')]);
  end
  disp(sprintf('---- dimension %d-D done ----', dim));
end
