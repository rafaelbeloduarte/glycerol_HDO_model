intervalo_t = 0:0.001:6;
% f_val = 0,3367
constante = [4.7461e-02
   2.5696e-02
   1.7791e-02
   4.4648e-02
   2.3021e-01
   4.8506e-01
   1.7238e-05];

  % Vetor Condições iniciais
  xG0 = 0.687; % mol/L
  Y0 = [0.00001 0.00001 0.00001 0.00001 0.00001 xG0 0.00001 transpose(constante)];

  t_exp = [0 1 2 3 4 5 6];
  x = [0.000120736019106005	0.00190542440690286	0.00126244275900532	0.00245650973286063	0.00345920725224455	0.0045534425575047	0.00547419411429379
0.000485247518719119	0.00231941595281414	0.00349664450705343	0.00495123526098734	0.00879527646171398	0.0128643049161132	0.0157702558358661
0.000313476087311038	0.00319273182508872	0.00613915891361497	0.0032510980857634	0.00474166685870984	0.00584445103395859	0.00737038552765253
0.686609344313597	0.674117033507055	0.646049321827957	0.619995278542439	0.578654902103528	0.566372569330732	0.534311242079215
0.000355852100700251	0.000452642255611177	0.00153275471395357	0.00313811173253462	0.007510876639749	0.011367323695359	0.017134530695938
0.00203921003920804	0.00549912899109952	0.0216181875618648	0.0372484791805704	0.054171111606807	0.0561776543209205	0.0675149706050606
0	0	0	0	0	0.0000102198652078027	0.0000176773328838207];
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

line_width = 2;
font_size = 16;
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
title("285^oC");
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
