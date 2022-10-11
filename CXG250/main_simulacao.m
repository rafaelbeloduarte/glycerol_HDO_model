intervalo_t = 0:0.001:6;
% f_val = 0,3367
constante = [1.8247e-02
   3.9749e-03
   6.5501e-03
   8.1070e-04
   2.3398e-02
   8.7818e-02
   8.4916e-04];

  % Vetor Condições iniciais
  xG0 = 0.683; % mol/L
  Y0 = [0.00001 0.00001 0.00001 0.00001 0.00001 xG0 0.00001 transpose(constante)];

  t_exp = [0 1 2 3 4 5 6];
  x = [0	0	0.00102032578230634	0.000847759680905678	0.00129073738355545	0.00154941781222605	0.00204699277401656
0.0000464469509095347	0.000145290674339949	0.000861137035061146	0.00231971482768176	0.00343554708446678	0.00451891816765057	0.00625638750654464
0	0	0.000158102767935209	0.000573856529033246	0.000842427332728351	0.000967574275187842	0.00131028807808756
0.682850238175542	0.681690995179061	0.677192383135557	0.670963636981612	0.664338057166123	0.659021461606686	0.657342101542777
0.000200714556852915	0.00018855465668809	0.000171049818282599	0	0.000209763576314985	0.000154137002728718	0.000250198727561956
0	0.000270784556533941	0.00146498626976766	0.00340171035416233	0.00463952695556679	0.00634848944545101	0.00703218851566364
0	0	0.000161043867173777	0.000715986773057563	0.000597574379113831	0.000922708411371734	0.00101437728916687];
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

    CI_lw = [1.71E-02
2.92E-03
5.44E-03
4.59E-08
2.22E-02
8.59E-02
7.13E-04
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
  CI_up = [2.02E-02
5.88E-03
8.42E-03
2.69E-03
2.53E-02
8.88E-02
1.34E-03
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

line_width = 2;
font_size = 14;
figure(1);
clf;
# 1,2-PDO
c = [0 0.443 0.737];
plot(t, x_12P, 'linewidth', line_width, 'color', c);
hold on;
scatter(t_exp, x_12P_exp, 20*line_width, c, "+", 'linewidth', line_width);
hold on;
plot(t, x_12P_lw, 'linewidth', line_width/line_width, 'linestyle', ":", 'color', c);
hold on;
plot(t, x_12P_up, 'linewidth', line_width/line_width, 'linestyle', ":", 'color', c);
hold on;
#1,3-PDO
c = [0.847 0.322 0.094];
plot(t, x_13P, 'linewidth', line_width, 'color', c);
hold on;
scatter(t_exp, x_13P_exp, 20*line_width, c, "o", 'linewidth', line_width);
hold on;
plot(t, x_13P_lw, "m", 'linewidth', line_width/line_width, 'linestyle', ":", 'color', c);
hold on;
plot(t, x_13P_up, "m", 'linewidth', line_width/line_width, 'linestyle', ":", 'color', c);
hold on;
# EG
c = [0.925 0.690 0.122];
plot(t, x_EG, 'linewidth', line_width, 'color', c);
hold on;
scatter(t_exp, x_EG_exp, 20*line_width, c, "*", 'linewidth', line_width);
hold on;
plot(t, x_EG_lw, "m", 'linewidth', line_width/line_width, 'linestyle', ":", 'color', c);
hold on;
plot(t, x_EG_up, "m", 'linewidth', line_width/line_width, 'linestyle', ":", 'color', c);
hold on;
# isobutil
c = [0.490 0.180 0.553];
plot(t, x_I, 'linewidth', line_width, 'color', c);
hold on;
scatter(t_exp, x_I_exp, 20*line_width, c, "^", 'linewidth', line_width);
hold on;
plot(t, x_I_lw, "m", 'linewidth', line_width/line_width, 'linestyle', ":", 'color', c);
hold on;
plot(t, x_I_up, "m", 'linewidth', line_width/line_width, 'linestyle', ":", 'color', c);
hold on;
# methane
c = [0.463 0.671 0.184];
plot(t, x_M, 'linewidth', line_width, 'color', c);
hold on;
scatter(t_exp, x_M_exp, 20*line_width, c, "v", 'linewidth', line_width);
hold on;
plot(t, x_M_lw, "m", 'linewidth', line_width/line_width, 'linestyle', ":", 'color', c);
hold on;
plot(t, x_M_up, "m", 'linewidth', line_width/line_width, 'linestyle', ":", 'color', c);
hold on;
# scopoletin
c = [0.298 0.741 0.929];
plot(t, x_S, 'linewidth', line_width, 'color', c);
hold on;
scatter(t_exp, x_S_exp, 20*line_width, c, "p", 'linewidth', line_width);
hold on;
plot(t, x_S_lw, "m", 'linewidth', line_width/line_width, 'linestyle', ":", 'color', c);
hold on;
plot(t, x_S_up, "m", 'linewidth', line_width/line_width, 'linestyle', ":", 'color', c);
hold on;

xlabel ("time (h)");
ylabel ("Molar fraction");
title("250^oC");
set(gca, "linewidth", line_width, "fontsize", font_size, 'fontweight', "bold", "fontname", "Liberation Serif", "ylim", [-inf inf]);
legend('x_{1,2-PDO}', 'x_{1,2-PDO, exp}',"x_{1,2-PDO, CI}","", 'x_{1,3-PDO}', 'x_{1,3-PDO, exp}',"x_{1,3-PDO, CI}","", 'x_{EG}', 'x_{EG,exp}',"x_{EG, CI}","", 'x_{I}', 'x_{I, exp}',"x_{I, CI}","", 'x_{CH_4}', 'x_{CH_4, exp}',"x_{CH_4, CI}","", 'x_{S}', 'x_{S, exp}',"x_{S, CI}","", 'fontsize', 14, 'fontweight', "bold", "linewidth", line_width,  'location', "eastoutside", "fontname", "Liberation Serif", "numcolumns", 1);
##set(gca, "linewidth", line_width, "fontsize", font_size, 'fontweight', "bold", "fontname", "Liberation Serif", "ylim", [-inf inf]);
##legend('x_{1,2-PDO}', 'x_{1,2-PDO, exp}',"","", 'x_{1,3-PDO}', 'x_{1,3-PDO, exp}',"","", 'x_{EG}', 'x_{EG,exp}',"","", 'x_{I}', 'x_{I, exp}',"","", 'x_{CH_4}', 'x_{CH_4, exp}',"","", 'x_{S}', 'x_{S, exp}',"","", 'fontsize', 24, 'fontweight', "bold", "linewidth", line_width,  'location', "eastoutside", "fontname", "Liberation Serif");
legend boxoff;
##linestyle
##
##    May be one of
##
##    ":"
##
##        Solid line. [default]
##    ":"
##
##        Dashed line.
##    ":"
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
plot(t, x_G_lw, 'k', 'linewidth', line_width/line_width, 'linestyle', ":");
hold on;
plot(t, x_G_up, 'k', 'linewidth', line_width/line_width, 'linestyle', ":");
hold on;
xlabel ("tempo (h)");
ylabel ("Fração molar");
set(gca, "linewidth", line_width, "fontsize", font_size);
legend('x_{Glicerol}', 'x_{Glicerol, exp}',"","", 'fontsize', font_size, 'location', "eastoutside");
