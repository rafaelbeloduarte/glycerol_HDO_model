% Função com o sistema de EDOs
function dX = odes_boot(t, X, k)
  k1 = k(1);
  k2 = k(2);
  k3 = k(3);
  k4 = k(4);
  k5 = k(5);
  k6 = k(6);
  k7 = k(7);
  V_L = 0.06; % L
  V_G = 0.3 - V_L; % L
  m_T = 1.26*V_L*1000; % g
  MM_G = 92.094;
  MM_12P = 76.095;
  MM_13P = 76.095;
  MM_EG = 62.068;
  MM_I = 158.24;
  MM_M = 16;
  MM_S = 192.16;
  % O vetor coluna das EDOs 'dX'
  % deve ser gerado pois a função 'edo45'aceita apenas
  % duas entradas no campo função (Variável independente, Variável dependente)
  % X(1) = P_CH3OH, X(2) = P_H2O, X(3) = P_H2, X(4) = P_CO2, X(5) = P_CO, X(6) = T
  % P_A = P_CH3OH, P_B = P_H2O, P_C = P_H2, P_D = P_CO2, P_E = P_CO, T=T0
  dX = zeros(7,1);
  n_T = m_T/(X(1)*MM_12P + X(2)*MM_13P + X(3)*MM_EG + X(4)*MM_S + X(5)*MM_I + X(6)*MM_G + X(7)*MM_M);
  dX(1) = V_L*r12P(k1)/n_T; % dx_12P/dt
  dX(2) = V_L*r13P(k2)/n_T; % dx_13P/dt
  dX(3) = V_L*rEG(k3)/n_T; % dx_EG/dt
  dX(4) = V_L*rS(k4)/n_T; % dx_S/dt
  dX(5) = V_L*rI(k5)/n_T; % dx_I/dt
  dX(6) = V_L*rG(k6)/n_T; % dx_G/dt
  dX(7) = V_G*rM(k7)/n_T; % dx_M/dt
end

function resultado_r12P = r12P(k1)
  r12P = k1;
  resultado_r12P = r12P;
end

function resultado_r13P = r13P(k2)
  r13P = k2;
  resultado_r13P = r13P;
end

function resultado_rEG = rEG(k3)
  rEG = k3;
  resultado_rEG = rEG;
end

function resultado_rS = rS(k4)
  rS = k4;
  resultado_rS = rS;
end

function resultado_rI = rI(k5)
  rI = k5;
  resultado_rI = rI;
end


function resultado_rG = rG(k6)
  rG = -k6;
  resultado_rG = rG;
end

function resultado_rM = rM(k7)
  rM = k7;
  resultado_rM = rM;
end
