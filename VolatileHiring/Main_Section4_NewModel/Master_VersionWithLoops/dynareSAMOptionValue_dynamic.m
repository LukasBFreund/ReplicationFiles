function [residual, g1, g2, g3] = dynareSAMOptionValue_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [M_.exo_nbr by nperiods] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(24, 1);
T19 = params(8)*y(7)^params(3);
T22 = y(8)^(1-params(3));
T162 = params(8)*getPowerDeriv(y(7),params(3),1);
T169 = getPowerDeriv(y(8),1-params(3),1);
lhs =y(7);
rhs =1-(1-params(4))*y(1);
residual(1)= lhs-rhs;
lhs =y(4);
rhs =T19*T22;
residual(2)= lhs-rhs;
lhs =y(11);
rhs =y(4)/y(7);
residual(3)= lhs-rhs;
lhs =y(9);
rhs =y(4)/y(8);
residual(4)= lhs-rhs;
lhs =y(5);
rhs =(1-params(4))*y(1)+y(4);
residual(5)= lhs-rhs;
lhs =y(6);
rhs =1-y(5);
residual(6)= lhs-rhs;
lhs =y(17);
rhs =(1-params(5))*params(9)*y(15)/(1-(1-params(4))*params(1))+y(21);
residual(7)= lhs-rhs;
lhs =y(21);
rhs =(1-params(5))*(params(9)*y(10)-params(7))+params(1)*(1-params(4))*y(29);
residual(8)= lhs-rhs;
lhs =y(18);
rhs =(1-params(5))*params(9)*y(16)/(1-(1-params(4))*params(1))+y(22);
residual(9)= lhs-rhs;
lhs =y(22);
rhs =(1-params(5))*(params(9)*y(10)-params(7))+params(1)*(1-params(4))*y(30);
residual(10)= lhs-rhs;
lhs =y(19);
rhs =params(1)*y(28)-y(14)*params(6)+y(9)*y(14)*(y(18)-params(1)*y(28));
residual(11)= lhs-rhs;
lhs =params(6);
rhs =y(9)*(y(17)-params(1)*y(28));
residual(12)= lhs-rhs;
lhs =y(14);
rhs =1-(y(15)-params(19))/(params(18)-params(19));
residual(13)= lhs-rhs;
lhs =y(16);
rhs =0.5*(y(15)+params(18));
residual(14)= lhs-rhs;
lhs =y(8);
rhs =y(14)*(params(16)-(1-params(4))*y(1));
residual(15)= lhs-rhs;
lhs =y(10);
rhs =(1-params(11))*(steady_state(7))+params(11)*y(2)+y(3)*x(it_, 1);
residual(16)= lhs-rhs;
lhs =y(13);
rhs =(1-params(12))*params(13)+y(3)*params(12)+params(14)*x(it_, 2);
residual(17)= lhs-rhs;
lhs =y(12);
rhs =y(5)*(y(10)+y(16));
residual(18)= lhs-rhs;
lhs =y(20);
rhs =(1-params(5))*params(9)*(steady_state(12))/(1-(1-params(4))*params(1))+y(24);
residual(19)= lhs-rhs;
lhs =y(24);
rhs =(1-params(5))*(params(9)*y(10)-params(7))+params(1)*(1-params(4))*y(31);
residual(20)= lhs-rhs;
lhs =y(23);
rhs =params(1)*y(28);
residual(21)= lhs-rhs;
lhs =y(25);
rhs =y(5);
residual(22)= lhs-rhs;
lhs =y(26);
rhs =y(10);
residual(23)= lhs-rhs;
lhs =y(27);
rhs =y(13);
residual(24)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(24, 33);

  %
  % Jacobian matrix
  %

  g1(1,1)=1-params(4);
  g1(1,7)=1;
  g1(2,4)=1;
  g1(2,7)=(-(T22*T162));
  g1(2,8)=(-(T19*T169));
  g1(3,4)=(-(1/y(7)));
  g1(3,7)=(-((-y(4))/(y(7)*y(7))));
  g1(3,11)=1;
  g1(4,4)=(-(1/y(8)));
  g1(4,8)=(-((-y(4))/(y(8)*y(8))));
  g1(4,9)=1;
  g1(5,4)=(-1);
  g1(5,1)=(-(1-params(4)));
  g1(5,5)=1;
  g1(6,5)=1;
  g1(6,6)=1;
  g1(7,15)=(-((1-params(5))*params(9)/(1-(1-params(4))*params(1))));
  g1(7,17)=1;
  g1(7,21)=(-1);
  g1(8,10)=(-((1-params(5))*params(9)));
  g1(8,21)=1;
  g1(8,29)=(-((1-params(4))*params(1)));
  g1(9,16)=(-((1-params(5))*params(9)/(1-(1-params(4))*params(1))));
  g1(9,18)=1;
  g1(9,22)=(-1);
  g1(10,10)=(-((1-params(5))*params(9)));
  g1(10,22)=1;
  g1(10,30)=(-((1-params(4))*params(1)));
  g1(11,9)=(-(y(14)*(y(18)-params(1)*y(28))));
  g1(11,14)=(-((-params(6))+y(9)*(y(18)-params(1)*y(28))));
  g1(11,18)=(-(y(9)*y(14)));
  g1(11,19)=1;
  g1(11,28)=(-(params(1)+y(9)*y(14)*(-params(1))));
  g1(12,9)=(-(y(17)-params(1)*y(28)));
  g1(12,17)=(-y(9));
  g1(12,28)=(-(y(9)*(-params(1))));
  g1(13,14)=1;
  g1(13,15)=1/(params(18)-params(19));
  g1(14,15)=(-0.5);
  g1(14,16)=1;
  g1(15,1)=(-(y(14)*(-(1-params(4)))));
  g1(15,8)=1;
  g1(15,14)=(-(params(16)-(1-params(4))*y(1)));
  g1(16,2)=(-params(11));
  g1(16,10)=1;
  g1(16,3)=(-x(it_, 1));
  g1(16,32)=(-y(3));
  g1(17,3)=(-params(12));
  g1(17,13)=1;
  g1(17,33)=(-params(14));
  g1(18,5)=(-(y(10)+y(16)));
  g1(18,10)=(-y(5));
  g1(18,12)=1;
  g1(18,16)=(-y(5));
  g1(19,20)=1;
  g1(19,24)=(-1);
  g1(20,10)=(-((1-params(5))*params(9)));
  g1(20,24)=1;
  g1(20,31)=(-((1-params(4))*params(1)));
  g1(21,28)=(-params(1));
  g1(21,23)=1;
  g1(22,5)=(-1);
  g1(22,25)=1;
  g1(23,10)=(-1);
  g1(23,26)=1;
  g1(24,13)=(-1);
  g1(24,27)=1;
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  v2 = zeros(32,3);
  v2(1,1)=2;
  v2(1,2)=205;
  v2(1,3)=(-(T22*params(8)*getPowerDeriv(y(7),params(3),2)));
  v2(2,1)=2;
  v2(2,2)=238;
  v2(2,3)=(-(T162*T169));
  v2(3,1)=2;
  v2(3,2)=206;
  v2(3,3)=  v2(2,3);
  v2(4,1)=2;
  v2(4,2)=239;
  v2(4,3)=(-(T19*getPowerDeriv(y(8),1-params(3),2)));
  v2(5,1)=3;
  v2(5,2)=202;
  v2(5,3)=(-((-1)/(y(7)*y(7))));
  v2(6,1)=3;
  v2(6,2)=106;
  v2(6,3)=  v2(5,3);
  v2(7,1)=3;
  v2(7,2)=205;
  v2(7,3)=(-((-((-y(4))*(y(7)+y(7))))/(y(7)*y(7)*y(7)*y(7))));
  v2(8,1)=4;
  v2(8,2)=235;
  v2(8,3)=(-((-1)/(y(8)*y(8))));
  v2(9,1)=4;
  v2(9,2)=107;
  v2(9,3)=  v2(8,3);
  v2(10,1)=4;
  v2(10,2)=239;
  v2(10,3)=(-((-((-y(4))*(y(8)+y(8))))/(y(8)*y(8)*y(8)*y(8))));
  v2(11,1)=11;
  v2(11,2)=438;
  v2(11,3)=(-(y(18)-params(1)*y(28)));
  v2(12,1)=11;
  v2(12,2)=278;
  v2(12,3)=  v2(11,3);
  v2(13,1)=11;
  v2(13,2)=570;
  v2(13,3)=(-y(14));
  v2(14,1)=11;
  v2(14,2)=282;
  v2(14,3)=  v2(13,3);
  v2(15,1)=11;
  v2(15,2)=575;
  v2(15,3)=(-y(9));
  v2(16,1)=11;
  v2(16,2)=447;
  v2(16,3)=  v2(15,3);
  v2(17,1)=11;
  v2(17,2)=900;
  v2(17,3)=(-(y(14)*(-params(1))));
  v2(18,1)=11;
  v2(18,2)=292;
  v2(18,3)=  v2(17,3);
  v2(19,1)=11;
  v2(19,2)=905;
  v2(19,3)=(-(y(9)*(-params(1))));
  v2(20,1)=11;
  v2(20,2)=457;
  v2(20,3)=  v2(19,3);
  v2(21,1)=12;
  v2(21,2)=537;
  v2(21,3)=(-1);
  v2(22,1)=12;
  v2(22,2)=281;
  v2(22,3)=  v2(21,3);
  v2(23,1)=12;
  v2(23,2)=900;
  v2(23,3)=params(1);
  v2(24,1)=12;
  v2(24,2)=292;
  v2(24,3)=  v2(23,3);
  v2(25,1)=15;
  v2(25,2)=430;
  v2(25,3)=1-params(4);
  v2(26,1)=15;
  v2(26,2)=14;
  v2(26,3)=  v2(25,3);
  v2(27,1)=16;
  v2(27,2)=1026;
  v2(27,3)=(-1);
  v2(28,1)=16;
  v2(28,2)=98;
  v2(28,3)=  v2(27,3);
  v2(29,1)=18;
  v2(29,2)=302;
  v2(29,3)=(-1);
  v2(30,1)=18;
  v2(30,2)=142;
  v2(30,3)=  v2(29,3);
  v2(31,1)=18;
  v2(31,2)=500;
  v2(31,3)=(-1);
  v2(32,1)=18;
  v2(32,2)=148;
  v2(32,3)=  v2(31,3);
  g2 = sparse(v2(:,1),v2(:,2),v2(:,3),24,1089);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  v3 = zeros(28,3);
  v3(1,1)=2;
  v3(1,2)=6739;
  v3(1,3)=(-(T22*params(8)*getPowerDeriv(y(7),params(3),3)));
  v3(2,1)=2;
  v3(2,2)=7828;
  v3(2,3)=(-(T169*params(8)*getPowerDeriv(y(7),params(3),2)));
  v3(3,1)=2;
  v3(3,2)=6740;
  v3(3,3)=  v3(2,3);
  v3(4,1)=2;
  v3(4,2)=6772;
  v3(4,3)=  v3(2,3);
  v3(5,1)=2;
  v3(5,2)=7861;
  v3(5,3)=(-(T162*getPowerDeriv(y(8),1-params(3),2)));
  v3(6,1)=2;
  v3(6,2)=6773;
  v3(6,3)=  v3(5,3);
  v3(7,1)=2;
  v3(7,2)=7829;
  v3(7,3)=  v3(5,3);
  v3(8,1)=2;
  v3(8,2)=7862;
  v3(8,3)=(-(T19*getPowerDeriv(y(8),1-params(3),3)));
  v3(9,1)=3;
  v3(9,2)=6736;
  v3(9,3)=(-((y(7)+y(7))/(y(7)*y(7)*y(7)*y(7))));
  v3(10,1)=3;
  v3(10,2)=3472;
  v3(10,3)=  v3(9,3);
  v3(11,1)=3;
  v3(11,2)=6640;
  v3(11,3)=  v3(9,3);
  v3(12,1)=3;
  v3(12,2)=6739;
  v3(12,3)=(-((y(7)*y(7)*y(7)*y(7)*(-(2*(-y(4))))-(-((-y(4))*(y(7)+y(7))))*(y(7)*y(7)*(y(7)+y(7))+y(7)*y(7)*(y(7)+y(7))))/(y(7)*y(7)*y(7)*y(7)*y(7)*y(7)*y(7)*y(7))));
  v3(13,1)=4;
  v3(13,2)=7858;
  v3(13,3)=(-((y(8)+y(8))/(y(8)*y(8)*y(8)*y(8))));
  v3(14,1)=4;
  v3(14,2)=3506;
  v3(14,3)=  v3(13,3);
  v3(15,1)=4;
  v3(15,2)=7730;
  v3(15,3)=  v3(13,3);
  v3(16,1)=4;
  v3(16,2)=7862;
  v3(16,3)=(-((y(8)*y(8)*y(8)*y(8)*(-(2*(-y(4))))-(-((-y(4))*(y(8)+y(8))))*(y(8)*y(8)*(y(8)+y(8))+y(8)*y(8)*(y(8)+y(8))))/(y(8)*y(8)*y(8)*y(8)*y(8)*y(8)*y(8)*y(8))));
  v3(17,1)=11;
  v3(17,2)=18951;
  v3(17,3)=(-1);
  v3(18,1)=11;
  v3(18,2)=9159;
  v3(18,3)=  v3(17,3);
  v3(19,1)=11;
  v3(19,2)=9287;
  v3(19,3)=  v3(17,3);
  v3(20,1)=11;
  v3(20,2)=14439;
  v3(20,3)=  v3(17,3);
  v3(21,1)=11;
  v3(21,2)=14727;
  v3(21,3)=  v3(17,3);
  v3(22,1)=11;
  v3(22,2)=18791;
  v3(22,3)=  v3(17,3);
  v3(23,1)=11;
  v3(23,2)=29841;
  v3(23,3)=params(1);
  v3(24,1)=11;
  v3(24,2)=9169;
  v3(24,3)=  v3(23,3);
  v3(25,1)=11;
  v3(25,2)=9617;
  v3(25,3)=  v3(23,3);
  v3(26,1)=11;
  v3(26,2)=14449;
  v3(26,3)=  v3(23,3);
  v3(27,1)=11;
  v3(27,2)=15057;
  v3(27,3)=  v3(23,3);
  v3(28,1)=11;
  v3(28,2)=29681;
  v3(28,3)=  v3(23,3);
  g3 = sparse(v3(:,1),v3(:,2),v3(:,3),24,35937);
end
end
