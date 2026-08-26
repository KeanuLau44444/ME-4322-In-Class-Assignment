% six bar linkage
% Static Equilibrium

clc;
clear;

% define the joints
A = [7 4 0];
B = [5 16 0];
C = [25 25 0];
D = [23 10 0];
E = [18 35 0];
F = [43  32 0];
G = [45 17 0];

% Define the lengths of the bars
lAB = norm(B - A);
lBC = norm(C - B);
lCD = norm(D - C);
lDE = norm(E - D);
lEF = norm(F - E);
lFG = norm(G - F);

% Weight of Each Link
WAB = [0 -1 0];
WBEC = [0 -1 0];
WCD = [0 -1 0];
WEF = [0 -1 0];
WFG = [0 -1 0];

% center of mass of each link
S1 = (A+B)/2;
S2 = (B+C+E)/3;
S3 = (C+D)/2;
S4 = (E+F)/2;
S5 = (F+G)/2;

syms FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin

% follow the pattern for the forces
ForceA = [FAx FAy 0]; 
ForceB = [FBx FBy 0];
ForceC = [FCx FCy 0];
ForceD = [FDx FDy 0];
ForceE = [FEx FEy 0];
ForceF = [FFx FFy 0];
ForceG = [FGx FGy 0];
InputTorque = [0 0 Tin];

% Applied Force
AppliedForce = [50 0 0];

% Static Equilibrium Condition for Link AB
eqn1 = ForceA + ForceB + WAB == 0;
eqn2 = cross(A-S1,ForceA) + cross(B-S1,ForceB) + InputTorque == 0;

eqn3 = -ForceB + ForceC + ForceE + WBEC == 0;
eqn4 = cross(B - S2, -ForceB) + cross(C - S2, ForceC) + cross(E - S2, ForceE) == 0;

eqn5 = -ForceC + ForceD + WCD == 0;
eqn6 = cross(C - S3, -ForceC) + cross(D - S3, ForceD) == 0;

eqn7 = -ForceE + ForceF + WEF == 0; 
eqn8 = cross(E - S4, -ForceE) + cross(F - S4, ForceF) == 0;

eqn9 = -ForceF + ForceG + WFG + AppliedForce == 0;
eqn10 = cross(F - S5, -ForceF) + cross(G - S5, ForceG) == 0;

% Solving the 10 equations
eqnMatrix = [eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8, eqn9, eqn10];

StaticSolution = solve(eqnMatrix, [FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin]);

Force_Ax = double(StaticSolution.FAx);
Force_Ay = double(StaticSolution.FAy);

Force_Bx = double(StaticSolution.FBx);
Force_By = double(StaticSolution.FBy);

Force_Cx = double(StaticSolution.FCx);
Force_Cy = double(StaticSolution.FCy);

Force_Dx = double(StaticSolution.FDx);
Force_Dy = double(StaticSolution.FDy);

Force_Ex = double(StaticSolution.FEx);
Force_Ey = double(StaticSolution.FEy);

Force_Fx = double(StaticSolution.FFx);
Force_Fy = double(StaticSolution.FFy);

Force_Gx = double(StaticSolution.FGx);
Force_Gy = double(StaticSolution.FGy);

Input_Torque = double(StaticSolution.Tin);


% Angular Velocity Calculation
% Loop ABCDA

syms wBEC wCD
omega_AB = [0 0 1];
omega_BEC = [0 0 wBEC];
omega_CD = [0 0 wCD];

eqn11 = cross(omega_AB,B-A) + cross(omega_BEC,C-B) + cross(omega_CD,D-C) == 0;

loopSolution = solve(eqn11,[wBEC wCD]);

angularVelocity_BEC = double(loopSolution.wBEC) %[output:7380f5ef]
angularVelocity_CD = double(loopSolution.wCD) %[output:72f9a00a]

omegaBEC = [0 0 angularVelocity_BEC];
omegaCD = [0 0 angularVelocity_CD];

syms wEF wFG

omega_EF = [0 0 wEF];
omega_GF = [0 0 wFG];

eqn12 = cross(omegaCD,C-D) + cross(omegaBEC,E-C) + cross(omega_EF,F-E) + cross(omega_GF,G-F) == 0;

loop2Solution = solve(eqn12, [wEF wFG]);

angularVelocity_EF = double(loop2Solution.wEF) %[output:330bc630]
angularVelocity_FG = double(loop2Solution.wFG) %[output:798ac654]

% Angular Acceleration
% Loop 1 ABCDA
syms aBEC aCD
alpha_AB = [0 0 0];
alpha_BEC = [0 0 aBEC];
alpha_CD = [0 0 aCD];

a_B_A = cross(alpha_AB,B-A) + cross(omega_AB, cross(omega_AB,B-A));
a_C_B = cross(alpha_BEC,C-B) + cross(omegaBEC, cross(omegaBEC,C-B));
a_D_C = cross(alpha_CD,D-C) + cross(omegaCD, cross(omegaCD,D-C));

eqn13 = a_B_A + a_C_B + a_D_C == 0;

loop1AccSolution = solve(eqn13,[aBEC aCD]);

alphaBEC = double(loop1AccSolution.aBEC) %[output:4c03174e]
alphaCD = double(loop1AccSolution.aCD) %[output:8e42b259]

alphaBEC_vector = [0 0 alphaBEC];
alphaCD_vector = [0 0 alphaCD];

syms aEF aFG
alpha_EF = [0 0 aEF];
alpha_FG = [0 0 aFG];

a_C_D = cross(alphaCD_vector,C-D) + cross(omegaCD,cross(omegaCD,C-D));
a_E_C = cross(alphaBEC_vector,E-C) + cross(omegaBEC,cross(omegaBEC,E-C));

angVEL_EF = [0 0 angularVelocity_EF];
angVEL_FG = [0 0 angularVelocity_FG];

a_E_F = cross(alpha_EF,F-E) + cross(angVEL_EF, cross(angVEL_EF,F-E));
a_G_F = cross(alpha_FG,G-F) + cross(angVEL_FG, cross(angVEL_FG,G-F));

eqn14 = a_C_D + a_E_C + a_E_F + a_G_F == 0;

loop2AccSolution = solve(eqn14, [aEF aFG]);

alphaEF = double(loop2AccSolution.aEF) %[output:788311f5]
alphaFG = double(loop2AccSolution.aFG) %[output:3fc79582]


% Velocity at Joint


% Velocity of S1

vB_A = cross(omega_AB,B-A);

% vS1_A = V_S1_B + V_B_A

V_S1_B = cross(omega_AB,S1-B);

vS1_A = V_S1_B + vB_A %[output:1750ede1]


% Velocity of S2

vB_A = cross(omega_AB,B-A);

% vS2_A = V_S2_B + V_B_A

V_S2_B = cross(omegaBEC,S2-B);

vS2_A = V_S2_B + vB_A %[output:7d75a9a1]


% Velocity of S3

vC_D = cross(omegaCD,C-D);

% vS3_D = V_S3_C + V_C_D

V_S3_C = cross(omegaCD,S3-C);

vS3_D = V_S3_C + vC_D %[output:44129a2a]


% Velocity of S4

V_F_G = cross(angVEL_FG,F-G);

% vS4_G = V_S4_F + V_F_G

V_S4_F = cross(angVEL_EF,S4-F);

vS4_G = V_S4_F + V_F_G %[output:0873e451]


% Velocity of S5

V_F_G = cross(angVEL_FG,F-G);

% vS5_G = V_S5_F + V_F_G

V_S5_F = cross(angVEL_FG,S5-F);

vS5_G = V_S5_F + V_F_G %[output:76369a13]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":23.2}
%---
%[output:7380f5ef]
%   data: {"dataType":"textualVariable","outputData":{"name":"angularVelocity_BEC","value":"0.1915"}}
%---
%[output:72f9a00a]
%   data: {"dataType":"textualVariable","outputData":{"name":"angularVelocity_CD","value":"0.9149"}}
%---
%[output:330bc630]
%   data: {"dataType":"textualVariable","outputData":{"name":"angularVelocity_EF","value":"-0.1047"}}
%---
%[output:798ac654]
%   data: {"dataType":"textualVariable","outputData":{"name":"angularVelocity_FG","value":"1.0635"}}
%---
%[output:4c03174e]
%   data: {"dataType":"textualVariable","outputData":{"name":"alphaBEC","value":"-0.0328"}}
%---
%[output:8e42b259]
%   data: {"dataType":"textualVariable","outputData":{"name":"alphaCD","value":"-0.2158"}}
%---
%[output:788311f5]
%   data: {"dataType":"textualVariable","outputData":{"name":"alphaEF","value":"-0.1596"}}
%---
%[output:3fc79582]
%   data: {"dataType":"textualVariable","outputData":{"name":"alphaFG","value":"0.0578"}}
%---
%[output:1750ede1]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"vS1_A","rows":1,"type":"double","value":[["-6","-1","0"]]}}
%---
%[output:7d75a9a1]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"vS2_A","rows":1,"type":"double","value":[["-13.7872","0.1064","0"]]}}
%---
%[output:44129a2a]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"vS3_D","rows":1,"type":"double","value":[["-6.8617","0.9149","0"]]}}
%---
%[output:0873e451]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"vS4_G","rows":1,"type":"double","value":[["-15.7953","-0.8188","0"]]}}
%---
%[output:76369a13]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"vS5_G","rows":1,"type":"double","value":[["-7.9761","-1.0635","0"]]}}
%---
