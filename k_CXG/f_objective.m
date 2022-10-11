function Y = f_objective(params,T,k)
  A = params(1);
  E = params(2);
  R = 8.314; % J/mol.K
  k_predicted = [];
  for i = 1:1:length(T)
    k_predicted(end+1) = A*exp(-E/(R*T(i)));
  endfor
  Y = sum(sum((k - k_predicted).^2));
endfunction
