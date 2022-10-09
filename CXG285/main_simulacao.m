intervalo_t = 0:0.001:6;
% f_val = 0,3367
constante = [1.5479e-01
   3.2578e-02
   3.6571e-02
   1.8304e-02
   1.7844e-01
   5.4667e-01
   4.3070e-03];

  % Vetor Condições iniciais
  xG0 = 0.689; % mol/L
  Y0 = [0.00001 0.00001 0.00001 0.00001 0.00001 xG0 0.00001 transpose(constante)];

  t_exp = [0 1 2 3 4 5 6];
  x = [0.00319405967940873	0.00164746168966925	0.00435408264677725	0.00589579091206648	0.00992218399149192	0.00937843951014123	0.00891449586718024
0.00115000920151936	0.00380478917780572	0.0124775300589467	0.0190923030315172	0.0433548830767874	0.0423689731542582	0.0394945343094663
0.00100554616685905	0.0014772555551388	0.0046016690585132	0.00506590213688732	0.00808372637745789	0.0078195303291812	0.00879400824960989
0.686869092252988	0.676493129296004	0.638162863496033	0.60625089458323	0.558821105508447	0.552794270849324	0.533776665933563
0.000341420115327105	0.000426603985522671	0.00108735858928218	0.00241640184470745	0.00347481766940102	0.00473875896015277	0.00651928920195604
0.00179131307316157	0.00831071357286639	0.0169422638012227	0.0308377822129207	0.0339635701481284	0.0439271905208685	0.0545461053884707
0	0.000174004664341329	0.000646946161425272	0.00201749537546178	0.00376287562664706	0.0051148183064752	0.00516357939171532];
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

    CI_lw = [1.47E-01
2.53E-02
2.90E-02
1.02E-02
1.71E-01
5.36E-01
3.54E-03
];
  Y0 = [0.00001 0.00001 0.00001 0.00001 0.00001 xG0 0.00001 transpose(CI_lw)];
  [t Y_lw] = ode45(@odes, intervalo_t, Y0);
  x_12P_lw = Y_lw(:,1);
  x_13P_lw = Y_lw(:,2);
  x_EG_lw = Y_lw(:,3);
  x_S_lw = Y_lw (:,4);
  x_I_lw = Y_lw(:,5);
  x_G_lw = Y_lw(:,6);
  x_M_lw = Y_lw(:,7);
  CI_up = [1.63E-01
4.14E-02
4.54E-02
2.70E-02
1.87E-01
5.55E-01
6.69E-03
];
  Y0 = [0.00001 0.00001 0.00001 0.00001 0.00001 xG0 0.00001 transpose(CI_up)];
  [t Y_up] = ode45(@odes, intervalo_t, Y0);
  x_12P_up = Y_up(:,1);
  x_13P_up = Y_up(:,2);
  x_EG_up = Y_up(:,3);
  x_S_up = Y_up (:,4);
  x_I_up = Y_up(:,5);
  x_G_up = Y_up(:,6);
  x_M_up = Y_up(:,7);

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
      if dif < 0.0001
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

line_width = 1;
font_size = 14;
figure(1);
clf;
# 1,2-PDO
plot(t, x_12P, "g", 'linewidth', line_width);
hold on;
scatter(t_exp, x_12P_exp, "g", "+", 'linewidth', line_width);
hold on;
plot(t, x_12P_lw, "g", 'linewidth', line_width/line_width, 'linestyle', ":");
hold on;
plot(t, x_12P_up, "g", 'linewidth', line_width/line_width, 'linestyle', "--");
hold on;
#1,3-PDO
plot(t, x_13P, "m", 'linewidth', line_width);
hold on;
scatter(t_exp, x_13P_exp, "m", "o", 'linewidth', line_width);
hold on;
plot(t, x_13P_lw, "m", 'linewidth', line_width/line_width, 'linestyle', ":");
hold on;
plot(t, x_13P_up, "m", 'linewidth', line_width/line_width, 'linestyle', "--");
hold on;
# EG
plot(t, x_EG, "b", 'linewidth', line_width);
hold on;
scatter(t_exp, x_EG_exp, "b", "*", 'linewidth', line_width);
hold on;
plot(t, x_EG_lw, "b", 'linewidth', line_width/line_width, 'linestyle', ":");
hold on;
plot(t, x_EG_up, "b", 'linewidth', line_width/line_width, 'linestyle', "--");
hold on;
# isobutil
plot(t, x_I, "r", 'linewidth', line_width);
hold on;
scatter(t_exp, x_I_exp, "r", "^", 'linewidth', line_width);
hold on;
plot(t, x_I_lw, "r", 'linewidth', line_width/line_width, 'linestyle', ":");
hold on;
plot(t, x_I_up, "r", 'linewidth', line_width/line_width, 'linestyle', "--");
hold on;
# methane
plot(t, x_M, "c", 'linewidth', line_width);
hold on;
scatter(t_exp, x_M_exp, "c", "v", 'linewidth', line_width);
hold on;
plot(t, x_M_lw, "c", 'linewidth', line_width/line_width, 'linestyle', ":");
hold on;
plot(t, x_M_up, "c", 'linewidth', line_width/line_width, 'linestyle', "--");
hold on;
# scopoletin
plot(t, x_S, "k",'linewidth', line_width);
hold on;
scatter(t_exp, x_S_exp, "k", "p", 'linewidth', line_width);
hold on;
plot(t, x_S_lw, "k", 'linewidth', line_width/line_width, 'linestyle', ":");
hold on;
plot(t, x_S_up, "k", 'linewidth', line_width/line_width, 'linestyle', "--");
hold on;

xlabel ("time (h)");
ylabel ("Molar fraction");
title("285^oC");
set(gca, "linewidth", line_width, "fontsize", font_size, 'fontweight', "bold", "fontname", "Liberation Serif", "ylim", [-inf inf]);
legend('x_{1,2-PDO}', 'x_{1,2-PDO, exp}',"x_{1,2-PDO, CI}","", 'x_{1,3-PDO}', 'x_{1,3-PDO, exp}',"x_{1,3-PDO, CI}","", 'x_{EG}', 'x_{EG,exp}',"x_{EG, CI}","", 'x_{I}', 'x_{I, exp}',"x_{I, CI}","", 'x_{CH_4}', 'x_{CH_4, exp}',"x_{CH_4, CI}","", 'x_{S}', 'x_{S, exp}',"x_{S, CI}","", 'fontsize', 10, 'fontweight', "bold", "linewidth", line_width,  'location', "eastoutside", "fontname", "Liberation Serif", "numcolumns", 2);
##set(gca, "linewidth", line_width, "fontsize", font_size, 'fontweight', "bold", "fontname", "Liberation Serif", "ylim", [-inf inf]);
##legend('x_{1,2-PDO}', 'x_{1,2-PDO, exp}',"","", 'x_{1,3-PDO}', 'x_{1,3-PDO, exp}',"","", 'x_{EG}', 'x_{EG,exp}',"","", 'x_{I}', 'x_{I, exp}',"","", 'x_{CH_4}', 'x_{CH_4, exp}',"","", 'x_{S}', 'x_{S, exp}',"","", 'fontsize', 24, 'fontweight', "bold", "linewidth", line_width,  'location', "eastoutside", "fontname", "Liberation Serif");
legend boxoff;
##linestyle
##
##    May be one of
##
##    "--"
##
##        Solid line. [default]
##    "--"
##
##        Dashed line.
##    "--"
##
##        Dotted line.
##    "-."
##
##        A dash-dot line.
##    "none"
##
##        No line. Points will still be marked using the current Marker Style.

figure(2);
clf;
plot(t, x_G, 'k', 'linewidth', line_width);
hold on;
scatter(t_exp, x_G_exp, 'k', 'linewidth', line_width);
hold on;
plot(t, x_G_lw, 'k', 'linewidth', line_width/line_width, 'linestyle', "--");
hold on;
plot(t, x_G_up, 'k', 'linewidth', line_width/line_width, 'linestyle', "--");
hold on;
xlabel ("tempo (h)");
ylabel ("Fração molar");
set(gca, "linewidth", line_width, "fontsize", font_size);
legend('x_{Glicerol}', 'x_{Glicerol, exp}',"","", 'fontsize', font_size, 'location', "eastoutside");
