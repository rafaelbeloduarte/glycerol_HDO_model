intervalo_t = 0:0.01:6;
% f_val = 0,3367
constante = [0.0711266099963456
0.0142443243265413
0.0193547518424003
0.00669429418027291
0.0951599505233869
0.283905370037338
0.0023480527582845
];

##  x_EG_exp = x(1,:);
##  x_12P_exp = x(2,:);
##  x_13P_exp = x(3,:);
##  x_G_exp = x(4,:);
##  x_S_exp = x (5,:);
##  x_I_exp = x(6,:);
##  x_M_exp = x(7,:);

##  dY(1) = V_L*r12P(Y)/n_T; % dx_12P/dt
##  dY(2) = V_L*r13P(Y)/n_T; % dx_13P/dt
##  dY(3) = V_L*rEG(Y)/n_T; % dx_EG/dt
##  dY(4) = V_L*rS(Y)/n_T; % dx_S/dt
##  dY(5) = V_L*rI(Y)/n_T; % dx_I/dt
##  dY(6) = V_L*rG(Y)/n_T; % dx_G/dt
##  dY(7) = V_G*rM(Y)/n_T; % dx_M/dt

  % Vetor Condições iniciais
  xG0 = 0.697; % mol/L
  Y0 = [0.00001 0.00001 0.00001 0.00001 0.00001 xG0 0.00001 transpose(constante)];

  t_exp = [0 1 2 3 4 5 6];
  x = [0.000854591538052104	0.000735741867627219	0.00194497996821644	0.00320424754168804	0.00396335808191175	0.00493929036908437	0.00569302635494354
0.000214868686104798	0.00131685656790968	0.00439722017907256	0.00941648668994965	0.0139112157053546	0.0175640855506391	0.0243574983801011
0	0.000465682969222304	0.00152324554357987	0.00226607982767692	0.002829008888551	0.00380209469987456	0.00409257132061907
0.697410535444331	0.693071469364175	0.674592916711882	0.652972973236454	0.638002912615113	0.62317263780838	0.609259693065981
0.000133923627881811	0.000258008988989206	0.000373434229557013	0.00078300660576858	0.00117820048555334	0.0017452268097832	0.00204860963506344
0.00125001129840931	0.00236740273605756	0.00760322502490318	0.0156877081170168	0.0207096967113278	0.0246597928580262	0.0278815585466593
0.000020913672637073	0.0000921069735930886	0.000374644140212069	0.00155680628740693	0.0015618261590607	0.00232242059337922	0.00334293743232316];
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

    CI_lw = [6.93E-02
1.24E-02
1.75E-02
4.85E-03
9.33E-02
2.77E-01
2.03E-03
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
  CI_up = [7.64E-02
1.97E-02
2.50E-02
1.23E-02
1.01E-01
2.85E-01
3.77E-03
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

exp_errors = [0.000629284550951304	5.86E-05	0.000129410796522941	5.72E-05	0.000579723930613297	0.00166788249898177	0.000124416745847602
0.000640037462557054	5.90E-05	0.000128053748053336	5.73E-05	0.000575870256522937	0.00168311665723299	0.00012376874797449
0.000627734092781036	5.93E-05	0.000128998783334368	5.74E-05	0.000577187760504662	0.00168926188008105	0.000128235570980332
0.000626713432564261	5.98E-05	0.000124579527185372	5.68E-05	0.000578565234986026	0.00165787380782516	0.000125384022217832
0.000642167829059896	5.91E-05	0.000128314060402038	5.62E-05	0.000586740480082407	0.00168949915153637	0.000125910182975607
0.000634335493621922	5.99E-05	0.000126986327897842	5.59E-05	0.000572187461397981	0.00167571011971594	0.000125242803157606
0.000644796619184721	5.86E-05	0.000129055536164625	5.65E-05	0.000566394235154414	0.00166668866204917	0.000127113540444958];

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
title("275^oC");
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
