intervalo_t = 0:0.001:6;
% f_val = 0,3367
constante = [1.5144e-03
   5.2215e-04
   2.8781e-04
   4.1290e-05
   2.5737e-02
   3.4119e-02
   0];

  % Vetor Condições iniciais
  xG0 = 0.684; % mol/L
  Y0 = [0.00001 0.00001 0.00001 0.00001 0.00001 xG0 0.00001 transpose(constante)];

  t_exp = [0 1 2 3 4 5 6];
  x = [0	0	0	0.0000705673656836223	0	0.000224334021026903	0.0000263653895395729
0	0.000049194610207348	0.00008038826006635	0.000167626216044625	0.00023443991334025	0.000444875774429263	0.000533540426327624
0	0	0	0	0	0.000230194856854255	0.000233529596418148
0.684442205351682	0.683663031866356	0.681776781563699	0.679726838611213	0.678165959149904	0.674577022915608	0.672741971574695
0	0	0	0	0	0.0000603036664514772	0
0	0.000700825416546395	0.0021718615139465	0.00346159006059108	0.00477848010634989	0.0065692398626041	0.00834782644138291
0	0	0	0	0	0	0];
  x_EG_exp = x(1,:);
  x_12P_exp = x(2,:);
  x_13P_exp = x(3,:);
  x_G_exp = x(4,:);
  x_S_exp = x (5,:);
  x_I_exp = x(6,:);
  x_M_exp = x(7,:);

  % Resolução das EDOs utilizando o método de Dormand-Prince de ordem 4
  [t Y] = ode45(@odes, intervalo_t, Y0);
  % Abrindo o vetor solução em suas componentes
  x_12P = [];
  x_13P = [];
  x_EG = [];
  x_S = [];
  x_I = [];
  x_G = [];
  x_M = [];

  x_12P = Y(:,1);
  x_13P = Y(:,2);
  x_EG = Y(:,3);
  x_S = Y (:,4);
  x_I = Y(:,5);
  x_G = Y(:,6);
  x_M = Y(:,7);

    res_x_12P = [];
  res_x_13P = [];
  res_x_EG = [];
  res_x_S = [];
  res_x_I = [];
  res_x_G = [];
  res_x_M = [];

  fit_x_12P = [];
  fit_x_13P = [];
  fit_x_EG = [];
  fit_x_S = [];
  fit_x_I = [];
  fit_x_G = [];
  fit_x_M = [];

  for i = 1:1:length(t_exp)
    for j = 1:1:length(intervalo_t)
      dif = abs(t_exp(i) - intervalo_t(j));
      if dif < 0.001
          res_x_12P(end+1) = x_12P(j) - x_12P_exp(i);
          fit_x_12P(end+1) = x_12P(j);
          res_x_13P(end+1) = x_13P(j) - x_13P_exp(i);
          fit_x_13P(end+1) = x_13P(j);
          res_x_EG(end+1) = x_EG(j) - x_EG_exp(i);
          fit_x_EG(end+1) = x_EG(j);
          res_x_S(end+1) = x_S(j) - x_S_exp(i);
          fit_x_S(end+1) = x_S(j);
          res_x_I(end+1) = x_I(j) - x_I_exp(i);
          fit_x_I(end+1) = x_I(j);
          res_x_G(end+1) = x_G(j) - x_G_exp(i);
          fit_x_G(end+1) = x_G(j);
          res_x_M(end+1) = x_M(j) - x_M_exp(i);
          fit_x_M(end+1) = x_M(j);
      endif
    endfor
endfor

  res_x_12P_ST = [];
  res_x_13P_ST = [];
  res_x_EG_ST = [];
  res_x_S_ST = [];
  res_x_I_ST = [];
  res_x_G_ST = [];
  res_x_M_ST = [];

  for i = 1:1:length(res_x_12P)
    res_x_12P_ST(end+1) = (res_x_12P(i)-mean(res_x_12P))/std(res_x_12P);
    res_x_13P_ST(end+1) = (res_x_13P(i)-mean(res_x_13P))/std(res_x_13P);
    res_x_EG_ST(end+1) = (res_x_EG(i)-mean(res_x_EG))/std(res_x_EG);
    res_x_S_ST(end+1) = (res_x_S(i)-mean(res_x_S))/std(res_x_S);
    res_x_I_ST(end+1) = (res_x_I(i)-mean(res_x_I))/std(res_x_I);
    res_x_G_ST(end+1) = (res_x_G(i)-mean(res_x_G))/std(res_x_G);
    res_x_M_ST(end+1) = (res_x_M(i)-mean(res_x_M))/std(res_x_M);
  endfor

residuals = [res_x_12P_ST, res_x_13P_ST, res_x_EG_ST, res_x_I_ST, res_x_M_ST, res_x_S_ST, res_x_G_ST];
fitted = [fit_x_12P, fit_x_13P, fit_x_EG, fit_x_I, fit_x_M, fit_x_S, fit_x_G];

csvwrite("residuals.csv", transpose(residuals));
csvwrite("fitted.csv", transpose(fitted));

figure(3);
clf;
scatter(fit_x_12P, res_x_12P, 'linewidth', line_width);
hold on;
scatter(fit_x_13P, res_x_13P, 'linewidth', line_width);
hold on;
scatter(fit_x_EG, res_x_EG, 'linewidth', line_width);
hold on;
scatter(fit_x_S, res_x_S, 'linewidth', line_width);
hold on;
scatter(fit_x_I, res_x_I, 'linewidth', line_width);
hold on;
#scatter(fit_x_G, res_x_G, 'linewidth', line_width);
#hold on;
scatter(fit_x_M, res_x_M, 'linewidth', line_width);
hold on;
xlabel ("Fitted molar fraction");
ylabel ("Residual molar fraction");
set(gca, "linewidth", line_width, "fontsize", font_size, "xaxislocation", "origin");
#legend('Residual_{1,2PDO}', 'fontsize', font_size, 'location', "eastoutside");


##marker
##    ‘+’	crosshair
##    ‘o’	circle
##    ‘*’	star
##    ‘.’	point
##    ‘x’	cross
##    ‘s’	square
##    ‘d’	diamond
##    ‘^’	upward-facing triangle
##    ‘v’	downward-facing triangle
##    ‘>’	right-facing triangle
##    ‘<’	left-facing triangle
##    ‘p’	pentagram
##    ‘h’	hexagram
% Gerando o gráfico

line_width = 2;
font_size = 16;
figure(1);
clf;
plot(t, x_12P, 'linewidth', line_width);
hold on;
scatter(t_exp, x_12P_exp, "+", 'linewidth', line_width);
hold on;
plot(t, x_13P, 'linewidth', line_width);
hold on;
scatter(t_exp, x_13P_exp, "o", 'linewidth', line_width);
hold on;
plot(t, x_EG, 'linewidth', line_width);
hold on;
scatter(t_exp, x_EG_exp, "*", 'linewidth', line_width);
hold on;
plot(t, x_I, 'linewidth', line_width);
hold on;
scatter(t_exp, x_I_exp, "^", 'linewidth', line_width);
hold on;
plot(t, x_M, 'linewidth', line_width);
hold on;
scatter(t_exp, x_M_exp, "v", 'linewidth', line_width);
hold on;
plot(t, x_S, 'linewidth', line_width);
hold on;
scatter(t_exp, x_S_exp, "p", 'linewidth', line_width);
hold on;
xlabel ("time (h)");
ylabel ("Molar fraction");
title("250^oC");
set(gca, "linewidth", line_width, "fontsize", font_size, 'fontweight', "bold", "fontname", "Liberation Serif");
legend('x_{1,2-PDO}', 'x_{1,2-PDO, exp}', 'x_{1,3-PDO}', 'x_{1,3-PDO, exp}', 'x_{EG}', 'x_{EG,exp}', 'x_{I}', 'x_{I, exp}', 'x_{CH_4}', 'x_{CH_4, exp}', 'x_{S}', 'x_{S, exp}', 'fontsize', font_size, 'fontweight', "bold", "linewidth", line_width,  'location', "eastoutside", "fontname", "Liberation Serif");
legend boxoff;

figure(2);
clf;
plot(t, x_G, 'linewidth', line_width);
hold on;
scatter(t_exp, x_G_exp, 'linewidth', line_width);
hold on;
xlabel ("tempo (h)");
ylabel ("Fração molar");
set(gca, "linewidth", line_width, "fontsize", font_size);
legend('x_{Glicerol}', 'x_{Glicerol, exp}', 'fontsize', font_size, 'location', "eastoutside");
