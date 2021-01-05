function [fval,info,exit_flag,DLIK,Hess,SteadyState,trend_coeff,Model,DynareOptions,BayesInfo,DynareResults] = dsge_likelihood(xparam1,DynareDataset,DatasetInfo,DynareOptions,Model,EstimatedParameters,BayesInfo,BoundsInfo,DynareResults,derivatives_info)
% Evaluates the posterior kernel of a dsge model using the specified
% kalman_algo; the resulting posterior includes the 2*pi constant of the
% likelihood function

%@info:
%! @deftypefn {Function File} {[@var{fval},@var{exit_flag},@var{ys},@var{trend_coeff},@var{info},@var{Model},@var{DynareOptions},@var{BayesInfo},@var{DynareResults},@var{DLIK},@var{AHess}] =} dsge_likelihood (@var{xparam1},@var{DynareDataset},@var{DynareOptions},@var{Model},@var{EstimatedParameters},@var{BayesInfo},@var{DynareResults},@var{derivatives_flag})
%! @anchor{dsge_likelihood}
%! @sp 1
%! Evaluates the posterior kernel of a dsge model.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item xparam1
%! Vector of doubles, current values for the estimated parameters.
%! @item DynareDataset
%! Matlab's structure describing the dataset (initialized by dynare, see @ref{dataset_}).
%! @item DynareOptions
%! Matlab's structure describing the options (initialized by dynare, see @ref{options_}).
%! @item Model
%! Matlab's structure describing the Model (initialized by dynare, see @ref{M_}).
%! @item EstimatedParamemeters
%! Matlab's structure describing the estimated_parameters (initialized by dynare, see @ref{estim_params_}).
%! @item BayesInfo
%! Matlab's structure describing the priors (initialized by dynare, see @ref{bayesopt_}).
%! @item DynareResults
%! Matlab's structure gathering the results (initialized by dynare, see @ref{oo_}).
%! @item derivates_flag
%! Integer scalar, flag for analytical derivatives of the likelihood.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item fval
%! Double scalar, value of (minus) the likelihood.
%! @item info
%! Double vector, second entry stores penalty, first entry the error code.
%! @table @ @code
%! @item info==0
%! No error.
%! @item info==1
%! The model doesn't determine the current variables uniquely.
%! @item info==2
%! MJDGGES returned an error code.
%! @item info==3
%! Blanchard & Kahn conditions are not satisfied: no stable equilibrium.
%! @item info==4
%! Blanchard & Kahn conditions are not satisfied: indeterminacy.
%! @item info==5
%! Blanchard & Kahn conditions are not satisfied: indeterminacy due to rank failure.
%! @item info==6
%! The jacobian evaluated at the deterministic steady state is complex.
%! @item info==19
%! The steadystate routine thrown an exception (inconsistent deep parameters).
%! @item info==20
%! Cannot find the steady state, info(4) contains the sum of square residuals (of the static equations).
%! @item info==21
%! The steady state is complex, info(4) contains the sum of square of imaginary parts of the steady state.
%! @item info==22
%! The steady has NaNs.
%! @item info==23
%! M_.params has been updated in the steadystate routine and has complex valued scalars.
%! @item info==24
%! M_.params has been updated in the steadystate routine and has some NaNs.
%! @item info==26
%! M_.params has been updated in the steadystate routine and has negative/0 values in loglinear model.
%! @item info==30
%! Ergodic variance can't be computed.
%! @item info==41
%! At least one parameter is violating a lower bound condition.
%! @item info==42
%! At least one parameter is violating an upper bound condition.
%! @item info==43
%! The covariance matrix of the structural innovations is not positive definite.
%! @item info==44
%! The covariance matrix of the measurement errors is not positive definite.
%! @item info==45
%! Likelihood is not a number (NaN).
%! @item info==46
%! Likelihood is a complex valued number.
%! @item info==47
%! Posterior kernel is not a number (logged prior density is NaN)
%! @item info==48
%! Posterior kernel is a complex valued number (logged prior density is complex).
%! @end table
%! @item exit_flag
%! Integer scalar, equal to zero if the routine return with a penalty (one otherwise).
%! @item DLIK
%! Vector of doubles, score of the likelihood.
%! @item AHess
%! Matrix of doubles, asymptotic hessian matrix.
%! @item SteadyState
%! Vector of doubles, steady state level for the endogenous variables.
%! @item trend_coeff
%! Matrix of doubles, coefficients of the deterministic trend in the measurement equation.
%! @item Model
%! Matlab's structure describing the model (initialized by dynare, see @ref{M_}).
%! @item DynareOptions
%! Matlab's structure describing the options (initialized by dynare, see @ref{options_}).
%! @item BayesInfo
%! Matlab's structure describing the priors (initialized by dynare, see @ref{bayesopt_}).
%! @item DynareResults
%! Matlab's structure gathering the results (initialized by dynare, see @ref{oo_}).
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{dynare_estimation_1}, @ref{mode_check}
%! @sp 2
%! @strong{This function calls:}
%! @sp 1
%! @ref{dynare_resolve}, @ref{lyapunov_symm}, @ref{lyapunov_solver}, @ref{compute_Pinf_Pstar}, @ref{kalman_filter_d}, @ref{missing_observations_kalman_filter_d}, @ref{univariate_kalman_filter_d}, @ref{kalman_steady_state}, @ref{getH}, @ref{kalman_filter}, @ref{score}, @ref{AHessian}, @ref{missing_observations_kalman_filter}, @ref{univariate_kalman_filter}, @ref{priordens}
%! @end deftypefn
%@eod:

% Copyright (C) 2004-2018 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

% AUTHOR(S) stephane DOT adjemian AT univ DASH lemans DOT FR

% Initialization of the returned variables and others...
fval        = [];
SteadyState = [];
trend_coeff = [];
exit_flag   = 1;
info        = zeros(4,1);
DLIK        = [];
Hess        = [];

% Ensure that xparam1 is a column vector.
xparam1 = xparam1(:);

% Set flag related to analytical derivatives.
analytic_derivation = DynareOptions.analytic_derivation;

if analytic_derivation && DynareOptions.loglinear
    error('The analytic_derivation and loglinear options are not compatible')
end

if nargout==1
    analytic_derivation=0;
end

if analytic_derivation
    kron_flag=DynareOptions.analytic_derivation_mode;
end

%------------------------------------------------------------------------------
% 1. Get the structural parameters & define penalties
%------------------------------------------------------------------------------

% Return, with endogenous penalty, if some parameters are smaller than the lower bound of the prior domain.
if isestimation(DynareOptions) && ~isequal(DynareOptions.mode_compute,1) && any(xparam1<BoundsInfo.lb)
    k = find(xparam1<BoundsInfo.lb);
    fval = Inf;
    exit_flag = 0;
    info(1) = 41;
    info(4)= sum((BoundsInfo.lb(k)-xparam1(k)).^2);
    if analytic_derivation
        DLIK=ones(length(xparam1),1);
    end
    return
end

% Return, with endogenous penalty, if some parameters are greater than the upper bound of the prior domain.
if isestimation(DynareOptions) && ~isequal(DynareOptions.mode_compute,1) && any(xparam1>BoundsInfo.ub)
    k = find(xparam1>BoundsInfo.ub);
    fval = Inf;
    exit_flag = 0;
    info(1) = 42;
    info(4)= sum((xparam1(k)-BoundsInfo.ub(k)).^2);
    if analytic_derivation
        DLIK=ones(length(xparam1),1);
    end
    return
end

% Get the diagonal elements of the covariance matrices for the structural innovations (Q) and the measurement error (H).
Model = set_all_parameters(xparam1,EstimatedParameters,Model);

Q = Model.Sigma_e;
H = Model.H;

% Test if Q is positive definite.
if ~issquare(Q) || EstimatedParameters.ncx || isfield(EstimatedParameters,'calibrated_covariances')
    [Q_is_positive_definite, penalty] = ispd(Q(EstimatedParameters.Sigma_e_entries_to_check_for_positive_definiteness,EstimatedParameters.Sigma_e_entries_to_check_for_positive_definiteness));
    if ~Q_is_positive_definite
        fval = Inf;
        exit_flag = 0;
        info(1) = 43;
        info(4) = penalty;
        return
    end
    if isfield(EstimatedParameters,'calibrated_covariances')
        correct_flag=check_consistency_covariances(Q);
        if ~correct_flag
            penalty = sum(Q(EstimatedParameters.calibrated_covariances.position).^2);
            fval = Inf;
            exit_flag = 0;
            info(1) = 71;
            info(4) = penalty;
            return
        end
    end
end

% Test if H is positive definite.
if ~issquare(H) || EstimatedParameters.ncn || isfield(EstimatedParameters,'calibrated_covariances_ME')
    [H_is_positive_definite, penalty] = ispd(H(EstimatedParameters.H_entries_to_check_for_positive_definiteness,EstimatedParameters.H_entries_to_check_for_positive_definiteness));
    if ~H_is_positive_definite
        fval = Inf;
        info(1) = 44;
        info(4) = penalty;
        exit_flag = 0;
        return
    end
    if isfield(EstimatedParameters,'calibrated_covariances_ME')
        correct_flag=check_consistency_covariances(H);
        if ~correct_flag
            penalty = sum(H(EstimatedParameters.calibrated_covariances_ME.position).^2);
            fval = Inf;
            exit_flag = 0;
            info(1) = 72;
            info(4) = penalty;
            return
        end
    end
end


%------------------------------------------------------------------------------
% 2. call model setup & reduction program
%------------------------------------------------------------------------------

% Linearize the model around the deterministic steady state and extract the matrices of the state equation (T and R).
[T,R,SteadyState,info,Model,DynareOptions,DynareResults] = dynare_resolve(Model,DynareOptions,DynareResults,'restrict');

% Return, with endogenous penalty when possible, if dynare_resolve issues an error code (defined in resol).
if info(1)
    if info(1) == 3 || info(1) == 4 || info(1) == 5 || info(1)==6 ||info(1) == 19 ||...
                info(1) == 20 || info(1) == 21 || info(1) == 23 || info(1) == 26 || ...
                info(1) == 81 || info(1) == 84 ||  info(1) == 85 ||  info(1) == 86
        %meaningful second entry of output that can be used
        fval = Inf;
        info(4) = info(2);
        exit_flag = 0;
        if analytic_derivation
            DLIK=ones(length(xparam1),1);
        end
        return
    else
        fval = Inf;
        info(4) = 0.1;
        exit_flag = 0;
        if analytic_derivation
            DLIK=ones(length(xparam1),1);
        end
        return
    end
end

% check endogenous prior restrictions
info=endogenous_prior_restrictions(T,R,Model,DynareOptions,DynareResults);
if info(1)
    fval = Inf;
    info(4)=info(2);
    exit_flag = 0;
    if analytic_derivation
        DLIK=ones(length(xparam1),1);
    end
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lnprior = priordens(xparam1,BayesInfo.pshape,BayesInfo.p6,BayesInfo.p7,BayesInfo.p3,BayesInfo.p4);

if isfield(DynareOptions,'irfs_match_estimation')
[fval,lik] = irf_matching(Model,DynareResults,DynareOptions,lnprior);
DLIK       = [-lnprior; lik(:)];
	   return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Define a vector of indices for the observed variables. Is this really usefull?...
BayesInfo.mf = BayesInfo.mf1;

% Define the constant vector of the measurement equation.
if DynareOptions.noconstant
    constant = zeros(DynareDataset.vobs,1);
else
    if DynareOptions.loglinear
        constant = log(SteadyState(BayesInfo.mfys));
    else
        constant = SteadyState(BayesInfo.mfys);
    end
end

% Define the deterministic linear trend of the measurement equation.
if BayesInfo.with_trend
    [trend_addition, trend_coeff]=compute_trend_coefficients(Model,DynareOptions,DynareDataset.vobs,DynareDataset.nobs);
    trend = repmat(constant,1,DynareDataset.nobs)+trend_addition;
else
    trend_coeff = zeros(DynareDataset.vobs,1);
    trend = repmat(constant,1,DynareDataset.nobs);
end

% Get needed informations for kalman filter routines.
start = DynareOptions.presample+1;
Z = BayesInfo.mf;           %selector for observed variables
no_missing_data_flag = ~DatasetInfo.missing.state;
mm = length(T);             %number of states
pp = DynareDataset.vobs;    %number of observables
rr = length(Q);             %number of shocks
kalman_tol = DynareOptions.kalman_tol;
diffuse_kalman_tol = DynareOptions.diffuse_kalman_tol;
riccati_tol = DynareOptions.riccati_tol;
Y = transpose(DynareDataset.data)-trend;

%------------------------------------------------------------------------------
% 3. Initial condition of the Kalman filter
%------------------------------------------------------------------------------
kalman_algo = DynareOptions.kalman_algo;


diffuse_periods = 0;
expanded_state_vector_for_univariate_filter=0;
singular_diffuse_filter = 0;
switch DynareOptions.lik_init
  case 1% Standard initialization with the steady state of the state equation.
    if kalman_algo~=2
        % Use standard kalman filter except if the univariate filter is explicitely choosen.
        kalman_algo = 1;
    end
    Pstar=lyapunov_solver(T,R,Q,DynareOptions);
    Pinf  = [];
    a     = zeros(mm,1);
    Zflag = 0;
  case 2% Initialization with large numbers on the diagonal of the covariance matrix if the states (for non stationary models).
    if kalman_algo ~= 2
        % Use standard kalman filter except if the univariate filter is explicitely choosen.
        kalman_algo = 1;
    end
    Pstar = DynareOptions.Harvey_scale_factor*eye(mm);
    Pinf  = [];
    a     = zeros(mm,1);
    Zflag = 0;
  case 3% Diffuse Kalman filter (Durbin and Koopman)
        % Use standard kalman filter except if the univariate filter is explicitely choosen.
    if kalman_algo == 0
        kalman_algo = 3;
    elseif ~((kalman_algo == 3) || (kalman_algo == 4))
        error(['The model requires Diffuse filter, but you specified a different Kalman filter. You must set options_.kalman_algo ' ...
               'to 0 (default), 3 or 4'])
    end
    [Pstar,Pinf] = compute_Pinf_Pstar(Z,T,R,Q,DynareOptions.qz_criterium);
    Z =zeros(length(BayesInfo.mf),size(T,1));
    for i = 1:length(BayesInfo.mf)
        Z(i,BayesInfo.mf(i))=1;
    end
    Zflag = 1;
    % Run diffuse kalman filter on first periods.
    if (kalman_algo==3)
        % Multivariate Diffuse Kalman Filter
        Pstar0 = Pstar; % store Pstar
        if no_missing_data_flag
            [dLIK,dlik,a,Pstar] = kalman_filter_d(Y, 1, size(Y,2), ...
                                                  zeros(mm,1), Pinf, Pstar, ...
                                                  kalman_tol, diffuse_kalman_tol, riccati_tol, DynareOptions.presample, ...
                                                  T,R,Q,H,Z,mm,pp,rr);
        else
            [dLIK,dlik,a,Pstar] = missing_observations_kalman_filter_d(DatasetInfo.missing.aindex,DatasetInfo.missing.number_of_observations,DatasetInfo.missing.no_more_missing_observations, ...
                                                              Y, 1, size(Y,2), ...
                                                              zeros(mm,1), Pinf, Pstar, ...
                                                              kalman_tol, diffuse_kalman_tol, riccati_tol, DynareOptions.presample, ...
                                                              T,R,Q,H,Z,mm,pp,rr);
        end
        diffuse_periods = length(dlik);
        if isinf(dLIK)
            % Go to univariate diffuse filter if singularity problem.
            singular_diffuse_filter = 1;
            Pstar = Pstar0;
        end
    end
    if singular_diffuse_filter || (kalman_algo==4)
        % Univariate Diffuse Kalman Filter
        if isequal(H,0)
            H1 = zeros(pp,1);
            mmm = mm;
        else
            if all(all(abs(H-diag(diag(H)))<1e-14))% ie, the covariance matrix is diagonal...
                H1 = diag(H);
                mmm = mm;
            else
                %Augment state vector (follows Section 6.4.3 of DK (2012))
                expanded_state_vector_for_univariate_filter=1;
                if Zflag
                    Z1=Z;
                else
                    Z1=zeros(pp,size(T,2));
                    for jz=1:length(Z)
                        Z1(jz,Z(jz))=1;
                    end
                end
                Z = [Z1, eye(pp)];
                Zflag=1;
                T = blkdiag(T,zeros(pp));
                Q = blkdiag(Q,H);
                R = blkdiag(R,eye(pp));
                Pstar = blkdiag(Pstar,H);
                Pinf  = blkdiag(Pinf,zeros(pp));
                H1 = zeros(pp,1);
                mmm   = mm+pp;

            end
        end

        [dLIK,dlik,a,Pstar] = univariate_kalman_filter_d(DatasetInfo.missing.aindex,...
                                                         DatasetInfo.missing.number_of_observations,...
                                                         DatasetInfo.missing.no_more_missing_observations, ...
                                                         Y, 1, size(Y,2), ...
                                                         zeros(mmm,1), Pinf, Pstar, ...
                                                         kalman_tol, diffuse_kalman_tol, riccati_tol, DynareOptions.presample, ...
                                                         T,R,Q,H1,Z,mmm,pp,rr);
        diffuse_periods = size(dlik,1);
    end
    if isnan(dLIK)
        fval = Inf;
        info(1) = 45;
        info(4) = 0.1;
        exit_flag = 0;
        return
    end

  case 4% Start from the solution of the Riccati equation.
    if kalman_algo ~= 2
        kalman_algo = 1;
    end
    if isequal(H,0)
        [err,Pstar] = kalman_steady_state(transpose(T),R*Q*transpose(R),transpose(build_selection_matrix(Z,mm,length(Z))));
    else
        [err,Pstar] = kalman_steady_state(transpose(T),R*Q*transpose(R),transpose(build_selection_matrix(Z,mm,length(Z))),H);
    end
    if err
        disp(['dsge_likelihood:: I am not able to solve the Riccati equation, so I switch to lik_init=1!']);
        DynareOptions.lik_init = 1;
        Pstar=lyapunov_solver(T,R,Q,DynareOptions);
    end
    Pinf  = [];
    a = zeros(mm,1);
    Zflag = 0;
  case 5            % Old diffuse Kalman filter only for the non stationary variables
    [eigenvect, eigenv] = eig(T);
    eigenv = diag(eigenv);
    nstable = length(find(abs(abs(eigenv)-1) > 1e-7));
    unstable = find(abs(abs(eigenv)-1) < 1e-7);
    V = eigenvect(:,unstable);
    indx_unstable = find(sum(abs(V),2)>1e-5);
    stable = find(sum(abs(V),2)<1e-5);
    nunit = length(eigenv) - nstable;
    Pstar = DynareOptions.Harvey_scale_factor*eye(nunit);
    if kalman_algo ~= 2
        kalman_algo = 1;
    end
    R_tmp = R(stable, :);
    T_tmp = T(stable,stable);
    Pstar_tmp=lyapunov_solver(T_tmp,R_tmp,Q,DynareOptions);
    Pstar(stable, stable) = Pstar_tmp;
    Pinf  = [];
    a = zeros(mm,1);
    Zflag = 0;
  otherwise
    error('dsge_likelihood:: Unknown initialization approach for the Kalman filter!')
end

if analytic_derivation
    offset = EstimatedParameters.nvx;
    offset = offset+EstimatedParameters.nvn;
    offset = offset+EstimatedParameters.ncx;
    offset = offset+EstimatedParameters.ncn;
    no_DLIK = 0;
    full_Hess = analytic_derivation==2;
    asy_Hess = analytic_derivation==-2;
    outer_product_gradient = analytic_derivation==-1;
    if asy_Hess
        analytic_derivation=1;
    end
    if outer_product_gradient
        analytic_derivation=1;
    end
    DLIK = [];
    AHess = [];
    iv = DynareResults.dr.restrict_var_list;
    if nargin<10 || isempty(derivatives_info)
        [A,B,nou,nou,Model,DynareOptions,DynareResults] = dynare_resolve(Model,DynareOptions,DynareResults);
        if ~isempty(EstimatedParameters.var_exo)
            indexo=EstimatedParameters.var_exo(:,1);
        else
            indexo=[];
        end
        if ~isempty(EstimatedParameters.param_vals)
            indparam=EstimatedParameters.param_vals(:,1);
        else
            indparam=[];
        end
        if full_Hess
            [dum, DT, DOm, DYss, dum2, D2T, D2Om, D2Yss] = getH(A, B, EstimatedParameters, Model,DynareResults,DynareOptions,kron_flag,indparam,indexo,iv);
            clear dum dum2;
        else
            [dum, DT, DOm, DYss] = getH(A, B, EstimatedParameters, Model,DynareResults,DynareOptions,kron_flag,indparam,indexo,iv);
        end
    else
        DT = derivatives_info.DT(iv,iv,:);
        DOm = derivatives_info.DOm(iv,iv,:);
        DYss = derivatives_info.DYss(iv,:);
        if isfield(derivatives_info,'full_Hess')
            full_Hess = derivatives_info.full_Hess;
        end
        if full_Hess
            D2T = derivatives_info.D2T;
            D2Om = derivatives_info.D2Om;
            D2Yss = derivatives_info.D2Yss;
        end
        if isfield(derivatives_info,'no_DLIK')
            no_DLIK = derivatives_info.no_DLIK;
        end
        clear('derivatives_info');
    end
    DYss = [zeros(size(DYss,1),offset) DYss];
    DH=zeros([length(H),length(H),length(xparam1)]);
    DQ=zeros([size(Q),length(xparam1)]);
    DP=zeros([size(T),length(xparam1)]);
    if full_Hess
        for j=1:size(D2Yss,1)
            tmp(j,:,:) = blkdiag(zeros(offset,offset), squeeze(D2Yss(j,:,:)));
        end
        D2Yss = tmp;
        D2H=sparse(size(D2Om,1),size(D2Om,2)); %zeros([size(H),length(xparam1),length(xparam1)]);
        D2P=sparse(size(D2Om,1),size(D2Om,2)); %zeros([size(T),length(xparam1),length(xparam1)]);
        jcount=0;
    end
    if DynareOptions.lik_init==1
        for i=1:EstimatedParameters.nvx
            k =EstimatedParameters.var_exo(i,1);
            DQ(k,k,i) = 2*sqrt(Q(k,k));
            dum =  lyapunov_symm(T,DOm(:,:,i),DynareOptions.lyapunov_fixed_point_tol,DynareOptions.qz_criterium,DynareOptions.lyapunov_complex_threshold,[],DynareOptions.debug);
            %         kk = find(abs(dum) < 1e-12);
            %         dum(kk) = 0;
            DP(:,:,i)=dum;
            if full_Hess
                for j=1:i
                    jcount=jcount+1;
                    dum =  lyapunov_symm(T,dyn_unvech(D2Om(:,jcount)),DynareOptions.lyapunov_fixed_point_tol,DynareOptions.qz_criterium,DynareOptions.lyapunov_complex_threshold,[],DynareOptions.debug);
                    %             kk = (abs(dum) < 1e-12);
                    %             dum(kk) = 0;
                    D2P(:,jcount)=dyn_vech(dum);
                    %             D2P(:,:,j,i)=dum;
                end
            end
        end
    end
    offset = EstimatedParameters.nvx;
    for i=1:EstimatedParameters.nvn
        k = EstimatedParameters.var_endo(i,1);
        DH(k,k,i+offset) = 2*sqrt(H(k,k));
        if full_Hess
            D2H(k,k,i+offset,i+offset) = 2;
        end
    end
    offset = offset + EstimatedParameters.nvn;
    if DynareOptions.lik_init==1
        for j=1:EstimatedParameters.np
            dum =  lyapunov_symm(T,DT(:,:,j+offset)*Pstar*T'+T*Pstar*DT(:,:,j+offset)'+DOm(:,:,j+offset),DynareOptions.lyapunov_fixed_point_tol,DynareOptions.qz_criterium,DynareOptions.lyapunov_complex_threshold,[],DynareOptions.debug);
            %         kk = find(abs(dum) < 1e-12);
            %         dum(kk) = 0;
            DP(:,:,j+offset)=dum;
            if full_Hess
                DTj = DT(:,:,j+offset);
                DPj = dum;
                for i=1:j+offset
                    jcount=jcount+1;
                    DTi = DT(:,:,i);
                    DPi = DP(:,:,i);
                    D2Tij = reshape(D2T(:,jcount),size(T));
                    D2Omij = dyn_unvech(D2Om(:,jcount));
                    tmp = D2Tij*Pstar*T' + T*Pstar*D2Tij' + DTi*DPj*T' + DTj*DPi*T' + T*DPj*DTi' + T*DPi*DTj' + DTi*Pstar*DTj' + DTj*Pstar*DTi' + D2Omij;
                    dum = lyapunov_symm(T,tmp,DynareOptions.lyapunov_fixed_point_tol,DynareOptions.qz_criterium,DynareOptions.lyapunov_complex_threshold,[],DynareOptions.debug);
                    %             dum(abs(dum)<1.e-12) = 0;
                    D2P(:,jcount) = dyn_vech(dum);
                    %             D2P(:,:,j+offset,i) = dum;
                end
            end
        end
    end
    if analytic_derivation==1
        analytic_deriv_info={analytic_derivation,DT,DYss,DOm,DH,DP,asy_Hess};
    else
        analytic_deriv_info={analytic_derivation,DT,DYss,DOm,DH,DP,D2T,D2Yss,D2Om,D2H,D2P};
        clear DT DYss DOm DP D2T D2Yss D2Om D2H D2P
    end
else
    analytic_deriv_info={0};
end

%------------------------------------------------------------------------------
% 4. Likelihood evaluation
%------------------------------------------------------------------------------

singularity_has_been_detected = false;
% First test multivariate filter if specified; potentially abort and use univariate filter instead
if ((kalman_algo==1) || (kalman_algo==3))% Multivariate Kalman Filter
    if no_missing_data_flag
        if DynareOptions.block
            [err, LIK] = block_kalman_filter(T,R,Q,H,Pstar,Y,start,Z,kalman_tol,riccati_tol, Model.nz_state_var, Model.n_diag, Model.nobs_non_statevar);
            mexErrCheck('block_kalman_filter', err);
        elseif DynareOptions.fast_kalman_filter
            if diffuse_periods
                %kalman_algo==3 requires no diffuse periods (stationary
                %observables) as otherwise FE matrix will not be positive
                %definite
                fval = Inf;
                info(1) = 55;
                info(4) = 0.1;
                exit_flag = 0;
                return
            end
            [LIK,lik] = kalman_filter_fast(Y,diffuse_periods+1,size(Y,2), ...
                                           a,Pstar, ...
                                           kalman_tol, riccati_tol, ...
                                           DynareOptions.presample, ...
                                           T,Q,R,H,Z,mm,pp,rr,Zflag,diffuse_periods, ...
                                           analytic_deriv_info{:});
        else
            [LIK,lik] = kalman_filter(Y,diffuse_periods+1,size(Y,2), ...
                                      a,Pstar, ...
                                      kalman_tol, riccati_tol, ...
                                      DynareOptions.rescale_prediction_error_covariance, ...
                                      DynareOptions.presample, ...
                                      T,Q,R,H,Z,mm,pp,rr,Zflag,diffuse_periods, ...
                                      analytic_deriv_info{:});
        end
    else
        if 0 %DynareOptions.block
            [err, LIK,lik] = block_kalman_filter(DatasetInfo.missing.aindex,DatasetInfo.missing.number_of_observations,DatasetInfo.missing.no_more_missing_observations,...
                                                 T,R,Q,H,Pstar,Y,start,Z,kalman_tol,riccati_tol, Model.nz_state_var, Model.n_diag, Model.nobs_non_statevar);
        else
            [LIK,lik] = missing_observations_kalman_filter(DatasetInfo.missing.aindex,DatasetInfo.missing.number_of_observations,DatasetInfo.missing.no_more_missing_observations,Y,diffuse_periods+1,size(Y,2), ...
                                                           a, Pstar, ...
                                                           kalman_tol, DynareOptions.riccati_tol, ...
                                                           DynareOptions.rescale_prediction_error_covariance, ...
                                                           DynareOptions.presample, ...
                                                           T,Q,R,H,Z,mm,pp,rr,Zflag,diffuse_periods);
        end
    end
    if analytic_derivation
        LIK1=LIK;
        LIK=LIK1{1};
        lik1=lik;
        lik=lik1{1};
    end
    if isinf(LIK)
        if DynareOptions.use_univariate_filters_if_singularity_is_detected
            singularity_has_been_detected = true;
            if kalman_algo == 1
                kalman_algo = 2;
            else
                kalman_algo = 4;
            end
        else
            fval = Inf;
            info(1) = 50;
            info(4) = 0.1;
            exit_flag = 0;
            return
        end
    else
        if DynareOptions.lik_init==3
            LIK = LIK + dLIK;
            if analytic_derivation==0 && nargout>3
                if ~singular_diffuse_filter
                    lik = [dlik; lik];
                else
                    lik = [sum(dlik,2); lik];
                end
            end
        end
    end
end

if (kalman_algo==2) || (kalman_algo==4)
    % Univariate Kalman Filter
    % resetting measurement error covariance matrix when necessary following DK (2012), Section 6.4.3                                                          %
    if isequal(H,0)
        H1 = zeros(pp,1);
        mmm = mm;
        if analytic_derivation
            DH = zeros(pp,length(xparam1));
        end
    else
        if all(all(abs(H-diag(diag(H)))<1e-14))% ie, the covariance matrix is diagonal...
            H1 = diag(H);
            mmm = mm;
            clear('tmp')
            if analytic_derivation
                for j=1:pp
                    tmp(j,:)=DH(j,j,:);
                end
                DH=tmp;
            end
        else
            if ~expanded_state_vector_for_univariate_filter
                Z1=zeros(pp,size(T,2));
                for jz=1:length(Z)
                    Z1(jz,Z(jz))=1;
                end
                Z = [Z1, eye(pp)];
                Zflag=1;
                T = blkdiag(T,zeros(pp));
                Q = blkdiag(Q,H);
                R = blkdiag(R,eye(pp));
                Pstar = blkdiag(Pstar,H);
                Pinf  = blkdiag(Pinf,zeros(pp));
                H1 = zeros(pp,1);
                Zflag=1;
            end
            mmm   = mm+pp;
            if singularity_has_been_detected
                a = zeros(mmm,1);
            elseif ~expanded_state_vector_for_univariate_filter
                a = [a; zeros(pp,1)];
            end
        end
    end
    if analytic_derivation
        analytic_deriv_info{5}=DH;
    end
    [LIK, lik] = univariate_kalman_filter(DatasetInfo.missing.aindex,DatasetInfo.missing.number_of_observations,DatasetInfo.missing.no_more_missing_observations,Y,diffuse_periods+1,size(Y,2), ...
                                          a,Pstar, ...
                                          DynareOptions.kalman_tol, ...
                                          DynareOptions.riccati_tol, ...
                                          DynareOptions.presample, ...
                                          T,Q,R,H1,Z,mmm,pp,rr,Zflag,diffuse_periods,analytic_deriv_info{:});
    if analytic_derivation
        LIK1=LIK;
        LIK=LIK1{1};
        lik1=lik;
        lik=lik1{1};
    end
    if DynareOptions.lik_init==3
        LIK = LIK+dLIK;
        if analytic_derivation==0 && nargout>3
            lik = [dlik; lik];
        end
    end
end

if analytic_derivation
    if no_DLIK==0
        DLIK = LIK1{2};
        %                 [DLIK] = score(T,R,Q,H,Pstar,Y,DT,DYss,DOm,DH,DP,start,Z,kalman_tol,riccati_tol);
    end
    if full_Hess
        Hess = -LIK1{3};
        %                     [Hess, DLL] = get_Hessian(T,R,Q,H,Pstar,Y,DT,DYss,DOm,DH,DP,D2T,D2Yss,D2Om,D2H,D2P,start,Z,kalman_tol,riccati_tol);
        %                     Hess0 = getHessian(Y,T,DT,D2T, R*Q*transpose(R),DOm,D2Om,Z,DYss,D2Yss);
    end
    if asy_Hess
        %         if ~((kalman_algo==2) || (kalman_algo==4)),
        %             [Hess] = AHessian(T,R,Q,H,Pstar,Y,DT,DYss,DOm,DH,DP,start,Z,kalman_tol,riccati_tol);
        %         else
        Hess = LIK1{3};
        %         end
    end
end

if isnan(LIK)
    fval = Inf;
    info(1) = 45;
    info(4) = 0.1;
    exit_flag = 0;
    return
end

if imag(LIK)~=0
    fval = Inf;
    info(1) = 46;
    info(4) = 0.1;
    exit_flag = 0;
    return
end

if isinf(LIK)~=0
    fval = Inf;
    info(1) = 50;
    info(4) = 0.1;
    exit_flag = 0;
    return
end

likelihood = LIK;

% ------------------------------------------------------------------------------
% 5. Adds prior if necessary
% ------------------------------------------------------------------------------
if analytic_derivation
    if full_Hess
        [lnprior, dlnprior, d2lnprior] = priordens(xparam1,BayesInfo.pshape,BayesInfo.p6,BayesInfo.p7,BayesInfo.p3,BayesInfo.p4);
        Hess = Hess - d2lnprior;
    else
        [lnprior, dlnprior] = priordens(xparam1,BayesInfo.pshape,BayesInfo.p6,BayesInfo.p7,BayesInfo.p3,BayesInfo.p4);
    end
    if no_DLIK==0
        DLIK = DLIK - dlnprior';
    end
    if outer_product_gradient
        dlik = lik1{2};
        dlik=[- dlnprior; dlik(start:end,:)];
        Hess = dlik'*dlik;
    end
else
    lnprior = priordens(xparam1,BayesInfo.pshape,BayesInfo.p6,BayesInfo.p7,BayesInfo.p3,BayesInfo.p4);
end

if DynareOptions.endogenous_prior==1
    if DynareOptions.lik_init==2 || DynareOptions.lik_init==3
        error('Endogenous prior not supported with non-stationary models')
    else
        [lnpriormom]  = endogenous_prior(Y,Pstar,BayesInfo,H);
        fval    = (likelihood-lnprior-lnpriormom);
    end
else
    fval    = (likelihood-lnprior);
end

if DynareOptions.prior_restrictions.status
    tmp = feval(DynareOptions.prior_restrictions.routine, Model, DynareResults, DynareOptions, DynareDataset, DatasetInfo);
    fval = fval - tmp;
end

if isnan(fval)
    fval = Inf;
    info(1) = 47;
    info(4) = 0.1;
    exit_flag = 0;
    return
end

if imag(fval)~=0
    fval = Inf;
    info(1) = 48;
    info(4) = 0.1;
    exit_flag = 0;
    return
end

if ~DynareOptions.kalman.keep_kalman_algo_if_singularity_is_detected
    % Update DynareOptions.kalman_algo.
    DynareOptions.kalman_algo = kalman_algo;
end

if analytic_derivation==0 && nargout>3
    lik=lik(start:end,:);
    DLIK=[-lnprior; lik(:)];
end
