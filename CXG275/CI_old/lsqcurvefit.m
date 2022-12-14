%{
This program is part of a set of scripts to fit a zeroth order reaction model
to experimental data as well as calculate the parameters confidence intervals.
Copyright (C) 2022  Rafael Belo Duarte

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

contact me at rafaelbeloduarte@pm.me
%}
pkg load optim;
 t_exp = [0 1 2 3 4 5 6];
  x_exp = transpose([0.000854591538052104	0.000735741867627219	0.00194497996821644	0.00320424754168804	0.00396335808191175	0.00493929036908437	0.00569302635494354
0.000214868686104798	0.00131685656790968	0.00439722017907256	0.00941648668994965	0.0139112157053546	0.0175640855506391	0.0243574983801011
0	0.000465682969222304	0.00152324554357987	0.00226607982767692	0.002829008888551	0.00380209469987456	0.00409257132061907
0.697410535444331	0.693071469364175	0.674592916711882	0.652972973236454	0.638002912615113	0.62317263780838	0.609259693065981
0.000133923627881811	0.000258008988989206	0.000373434229557013	0.00078300660576858	0.00117820048555334	0.0017452268097832	0.00204860963506344
0.00125001129840931	0.00236740273605756	0.00760322502490318	0.0156877081170168	0.0207096967113278	0.0246597928580262	0.0278815585466593
0.000020913672637073	0.0000921069735930886	0.000374644140212069	0.00155680628740693	0.0015618261590607	0.00232242059337922	0.00334293743232316]);

k0 = [1
   4.2321e-02
   1
   6.9700e-01
   1
   3.0187e-02
   6.8913e-04];
##
##%% bounds
##lb=[];
##ub=[];
##
##[k, resnorm, residual, flag, output, lambda, jacob] = lsqnonlin(@(x) odes_nova(t_exp, x, x_exp, x0), k0, lb,  ub)

xG0 = 0.697; % mol/L
init = [0.00001 0.00001 0.00001 xG0 0.00001 0.00001 0.00001];
newparam=lsqnonlin(@(k)fit_simp(k,t_exp,x_exp),k0,init,[])
##  x_EG_exp = x(1,:);
##  x_12P_exp = x(2,:);
##  x_13P_exp = x(3,:);
##  x_G_exp = x(4,:);
##  x_S_exp = x (5,:);
##  x_I_exp = x(6,:);
##  x_M_exp = x(7,:);

