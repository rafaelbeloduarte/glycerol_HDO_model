intervalo_t = 0:0.001:6;
% f_val = 0,3367
constante = [1.1303e-02
   5.8348e-03
   5.9348e-03
   3.5684e-03
   8.6974e-02
   1.5845e-01
   1.4227e-05];

  % Vetor Condições iniciais
  xG0 = 0.677; % mol/L
  Y0 = [0.00001 0.00001 0.00001 0.00001 0.00001 xG0 0.00001 transpose(constante)];

  t_exp = [0 1 2 3 4 5 6];
  x = [0	0.000258008137456242	0.00121734948945218	0.000990533154712822	0.00151925377890684	0.00189741662954417	0.000962807544541986
0	0.000347947079806869	0.00117345315677939	0.00171019875794062	0.00272965691238154	0.00293876797155187	0.00296712067296808
0	0.000349178250508359	0.000884365384512156	0.000918912902739483	0.0016725747998978	0.00118313576345039	0.00151098201015375
0.677031024434055	0.670844940635314	0.663018117833216	0.649846337997237	0.643460785061871	0.642518115451164	0.628074600961928
0	0.000238936100932408	0.000254443701267073	0.000502359904513105	0.000625801411513753	0.000806059076793232	0.00124700350449904
0	0.0033193723371167	0.00846468693777631	0.012777414103993	0.0165508380532851	0.0213796891746793	0.0268823418573825
0	0	0	0	0.0000238830013814539	0.0000291758702947773	0.0000369245562455555];
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
title("265^oC");
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
