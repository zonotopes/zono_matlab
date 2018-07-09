%% Example of zonotope library call for a dataset
%% of 10 generators in 3 dimensions

star = [ 9.7974838e-01   6.2406009e-01   6.0986665e-01 ;
         4.3886997e-01   6.7913554e-01   6.1766639e-01 ;
         1.1111922e-01   3.9551522e-01   8.5944231e-01 ;
         2.5806470e-01   3.6743665e-01   8.0548942e-01 ;
         4.0871985e-01   9.8798200e-01   5.7672152e-01 ;
         5.9489607e-01   3.7738866e-02   1.8292247e-01 ;
         2.6221175e-01   8.8516801e-01   2.3993201e-01 ;
         6.0284309e-01   9.1328683e-01   8.8651193e-01 ;
         7.1121578e-01   7.9618387e-01   2.8674152e-02 ;
         2.2174673e-01   9.8712279e-02   4.8990139e-01 ];
     
computeAdditionalStatistics = true; % put this flag to false if your
                                    % are interested in volume only
                                    
volume = zonotope(star, computeAdditionalStatistics);

fprintf('\nThe computed volume is %g\n', volume);
fprintf('(the known volume for this dataset is 20.1948)\n');