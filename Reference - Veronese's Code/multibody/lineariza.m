function A = lineariza(fun, vec)
  Ncol = length(vec);
  Nlin = length(fun(vec));
  delta = 1e-6;
  A = zeros(Nlin, Ncol);
  for  i = 1:Ncol
      Delta = zeros(Ncol,1);
      Delta(i) = delta;
      A(:,i) = (fun(vec+Delta) - fun(vec-Delta))/(2*delta);
  end
end