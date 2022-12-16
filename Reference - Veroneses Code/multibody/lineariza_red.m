function [A] = lineariza_red(t,eta0,n,Kn)
Ncol = length(eta0);
Nlin = length(dinam_aero_red(t,eta0,n,Kn));
delta = 1e-6;
A = zeros(Nlin, Ncol);
  for  i = 1:Ncol
      Delta = zeros(Ncol,1);
      Delta(i) = delta;
      A(:,i) = (dinam_aero_red(t,eta0+Delta,n,Kn) - dinam_aero_red(t,eta0-Delta,n,Kn))/(2*delta);
  end
end