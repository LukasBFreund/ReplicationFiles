function [residual, g1, g2] = dynareSAMOptionValue_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                     columns: variables in declaration order
%                                                     rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 24, 1);

%
% Model equations
%

T183 = (-((1-params(5))*params(9)/(1-(1-params(4))*params(1))));
lhs =y(4);
rhs =1-(1-params(4))*y(2);
residual(1)= lhs-rhs;
lhs =y(1);
rhs =params(8)*y(4)^params(3)*y(5)^(1-params(3));
residual(2)= lhs-rhs;
lhs =y(8);
rhs =y(1)/y(4);
residual(3)= lhs-rhs;
lhs =y(6);
rhs =y(1)/y(5);
residual(4)= lhs-rhs;
lhs =y(2);
rhs =(1-params(4))*y(2)+y(1);
residual(5)= lhs-rhs;
lhs =y(3);
rhs =1-y(2);
residual(6)= lhs-rhs;
lhs =y(14);
rhs =(1-params(5))*params(9)*y(12)/(1-(1-params(4))*params(1))+y(15);
residual(7)= lhs-rhs;
lhs =y(15);
rhs =(1-params(5))*(params(9)*y(7)-params(7))+params(1)*(1-params(4))*y(15);
residual(8)= lhs-rhs;
lhs =y(16);
rhs =(1-params(5))*params(9)*y(13)/(1-(1-params(4))*params(1))+y(17);
residual(9)= lhs-rhs;
lhs =y(17);
rhs =(1-params(5))*(params(9)*y(7)-params(7))+params(1)*(1-params(4))*y(17);
residual(10)= lhs-rhs;
lhs =y(18);
rhs =params(1)*y(18)+y(11)*((-params(6))/y(6)+y(16)-params(1)*y(18));
residual(11)= lhs-rhs;
lhs =params(6);
rhs =y(6)*(y(14)-params(1)*y(18));
residual(12)= lhs-rhs;
lhs =y(11);
rhs =1-(y(12)-params(19))/(params(18)-params(19));
residual(13)= lhs-rhs;
lhs =y(13);
rhs =0.5*(y(12)+params(18));
residual(14)= lhs-rhs;
lhs =y(5);
rhs =y(11)*(params(16)-(1-params(4))*y(2));
residual(15)= lhs-rhs;
lhs =y(7);
rhs =(1-params(11))*(y(7))+y(7)*params(11)+y(10)*x(1);
residual(16)= lhs-rhs;
lhs =y(10);
rhs =(1-params(12))*params(13)+y(10)*params(12)+params(14)*x(2);
residual(17)= lhs-rhs;
lhs =y(9);
rhs =y(2)*(y(7)+y(13));
residual(18)= lhs-rhs;
lhs =y(20);
rhs =(1-params(5))*params(9)*(y(12))/(1-(1-params(4))*params(1))+y(21);
residual(19)= lhs-rhs;
lhs =y(21);
rhs =(1-params(5))*(params(9)*y(7)-params(7))+params(1)*(1-params(4))*y(21);
residual(20)= lhs-rhs;
lhs =y(19);
rhs =params(1)*y(18);
residual(21)= lhs-rhs;
lhs =y(22);
rhs =y(2);
residual(22)= lhs-rhs;
lhs =y(23);
rhs =y(7);
residual(23)= lhs-rhs;
lhs =y(24);
rhs =y(10);
residual(24)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(24, 24);

  %
  % Jacobian matrix
  %

  g1(1,2)=1-params(4);
  g1(1,4)=1;
  g1(2,1)=1;
  g1(2,4)=(-(y(5)^(1-params(3))*params(8)*getPowerDeriv(y(4),params(3),1)));
  g1(2,5)=(-(params(8)*y(4)^params(3)*getPowerDeriv(y(5),1-params(3),1)));
  g1(3,1)=(-(1/y(4)));
  g1(3,4)=(-((-y(1))/(y(4)*y(4))));
  g1(3,8)=1;
  g1(4,1)=(-(1/y(5)));
  g1(4,5)=(-((-y(1))/(y(5)*y(5))));
  g1(4,6)=1;
  g1(5,1)=(-1);
  g1(5,2)=1-(1-params(4));
  g1(6,2)=1;
  g1(6,3)=1;
  g1(7,12)=T183;
  g1(7,14)=1;
  g1(7,15)=(-1);
  g1(8,7)=(-((1-params(5))*params(9)));
  g1(8,15)=1-(1-params(4))*params(1);
  g1(9,13)=T183;
  g1(9,16)=1;
  g1(9,17)=(-1);
  g1(10,7)=(-((1-params(5))*params(9)));
  g1(10,17)=1-(1-params(4))*params(1);
  g1(11,6)=(-(y(11)*params(6)/(y(6)*y(6))));
  g1(11,11)=(-((-params(6))/y(6)+y(16)-params(1)*y(18)));
  g1(11,16)=(-y(11));
  g1(11,18)=1-(params(1)+y(11)*(-params(1)));
  g1(12,6)=(-(y(14)-params(1)*y(18)));
  g1(12,14)=(-y(6));
  g1(12,18)=(-(y(6)*(-params(1))));
  g1(13,11)=1;
  g1(13,12)=1/(params(18)-params(19));
  g1(14,12)=(-0.5);
  g1(14,13)=1;
  g1(15,2)=(-(y(11)*(-(1-params(4)))));
  g1(15,5)=1;
  g1(15,11)=(-(params(16)-(1-params(4))*y(2)));
  g1(16,7)=1-(params(11)+1-params(11));
  g1(16,10)=(-x(1));
  g1(17,10)=1-params(12);
  g1(18,2)=(-(y(7)+y(13)));
  g1(18,7)=(-y(2));
  g1(18,9)=1;
  g1(18,13)=(-y(2));
  g1(19,12)=T183;
  g1(19,20)=1;
  g1(19,21)=(-1);
  g1(20,7)=(-((1-params(5))*params(9)));
  g1(20,21)=1-(1-params(4))*params(1);
  g1(21,18)=(-params(1));
  g1(21,19)=1;
  g1(22,2)=(-1);
  g1(22,22)=1;
  g1(23,7)=(-1);
  g1(23,23)=1;
  g1(24,10)=(-1);
  g1(24,24)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],24,576);
end
end
