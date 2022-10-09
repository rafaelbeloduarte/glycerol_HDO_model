intervalo_t = 0:0.001:6;
% f_val = 0,3367
constante = [2.9630e-02
   1.5402e-02
   1.5835e-02
   1.4668e-02
   1.6155e-01
   2.5903e-01
   9.0510e-06];

  % Vetor Condições iniciais
  xG0 = 0.689; % mol/L
  Y0 = [0.00001 0.00001 0.00001 0.00001 0.00001 xG0 0.00001 transpose(constante)];

  t_exp = [0 1 2 3 4 5 6];
  x = [0.00194549439845336	0.00087395975525803	0.00149818889892486	0.0024908991976411	0.00504508604594008	0.00387379954231023	0.00362180474518346
0.00216965789727529	0.00137752474610415	0.00352074784393782	0.00443059632909249	0.00664672827892196	0.00707658950004113	0.00881691751137682
0.00120648035047698	0.000741429597178869	0.00258519434767493	0.00311480245840776	0.00424196114852012	0.00333162652721618	0.00356864835517883
0.688822626364104	0.695865769598272	0.670716222663411	0.656100405778356	0.628654017166774	0.624994356029778	0.602407417663549
0	0.0000916192956785067	0.000626363388483587	0.000975518345421958	0.00223768655071472	0.00307931528475682	0.00588481057466947
0.000637319447767249	0.00461144893058251	0.0125505774448675	0.019782580977562	0.0324970220260622	0.0409478188712749	0.0542636268256832
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
title("275^oC");
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
