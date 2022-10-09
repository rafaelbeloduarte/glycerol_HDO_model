function Y = f_diff_sqr_boot(k, X_calculated, t_range)
  % Vetor Condições iniciais
  xG0 = 0.663; % mol/L
  X0 = [0.00001, 0.00001, 0.00001, 0.00001, 0.00001, xG0, 0.00001];
  % Resolução das EDOs utilizando o método de Dormand-Prince de ordem 4
  [t X] = ode45(@(t, X) odes_boot(t,X,abs(k)), t_range, X0);
  % Abrindo o vetor solução em suas componentes
  Y = sum(sum((X - X_calculated).^2));
endfunction

% [xmin, fval]=fminsearch(@f_250,[11.3], optimset ("Display", "iter"))
% [x, fval, info, output, grad, hess] = fminunc(@f_300,[10;0.5;0.1], optimset ("TolFun",0.0001))
