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

    CI_lw = [0.0092650036451264
0.00398275618117765
0.00387553211883219
0.0014286973792167
0.0851243216965927
0.156333069269707
2.28E-10
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
  CI_up = [0.0133884667429858
0.00794331183324324
0.00808761242765737
0.00567558351484677
0.0889857385974245
0.1604830008869
0.000653532911400163
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
title("265^oC");
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
