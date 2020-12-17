%==========================================================================
% Parent file: SaM model with option value due to het. prod. draws.
%
% This file: compute steady-state elasticity of hiring rate wrt agg prod, 
% and does so *numerically* (used for normal dist., here as a validity check)
%
%
% Case where y = z+a and entrepreneurs have 0 value upon separation
%
% Structure: 
%           Inputs: 
%                   - x: 2x1 vector of omega and chi guesses
%                   - sPar:   structure of parameter values
%
%           Functions called:
%                  - fn_Cutoff                
%
%           Output:  
%                   - elasticity: 1x1 scalar elasticity
%
% Last updated: November 2020
%==========================================================================

function elasticity = fn_Elasticity_Numerical(x,sPar)
    
    sParGuess = sPar;

    % Impose conjectured values of omega and chi 
    sParGuess.omega = x(1);
    sParGuess.chi = x(2);

    % Set up *numerical* computation of elasticity (here done for a 1% increase
    % in z; robustness and detailed results available)
    difz     = 0.001;
    dlnz     = log((sParGuess.ss.z+difz)/sParGuess.ss.z);
    sParGuess.ss.z = sParGuess.ss.z + difz;
  
    % Compute h after change in z
    vInit = [0.001,0.7];  
    options = optimoptions('fsolve','Display','none',...
                'MaxFunctionEvaluations',100000,'StepTolerance',1e-10,'OptimalityTolerance',1e-10);
    
     vSol = fsolve(@(x) fn_Cutoff(x,sParGuess),vInit,options);
     vSol = real(vSol);
     h = vSol(2);
    
    % Implied elasticity         
    elasticity =     log(h/sParGuess.ss.h)/dlnz;
   
end
