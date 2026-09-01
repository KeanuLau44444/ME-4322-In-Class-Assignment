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
F = [43 32 0];
G = [45 17 0];

% Joint values
new_B_x(1) = B(1);
new_B_y(1) = B(1);
new_E_x(1) = E(1);
new_E_y(1) = E(1);
new_F_x(1) = F(1);
new_F_y(1) = F(1);
   

% Define the lengths of the bars

lAB = norm(B - A);
lBC = norm(C - B);
lCD = norm(D - C);
lDE = norm(E - D);
lEF = norm(F - E);
lFG = norm(G - F);
lCE = norm(E - C);
lBE = norm(E - B);

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

eqn4 = cross(B-S2,-ForceB) + cross(C-S2,ForceC) + cross(E-S2,ForceE) == 0;

eqn5 = -ForceC + ForceD + WCD == 0;

eqn6 = cross(C-S3,-ForceC) + cross(D-S3,ForceD) == 0;

eqn7 = -ForceE + ForceF + WEF == 0;

eqn8 = cross(E-S4,-ForceE) + cross(F-S4,ForceF) == 0;

eqn9 = -ForceF + ForceG + WFG + AppliedForce == 0;

eqn10 = cross(F-S5,-ForceF) + cross(G-S5,ForceG) == 0;

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


%% Angular Velocity Calculation

% Loop 1 ABCDA

syms wBEC wCD

omega_AB = [0 0 1];

omega_BEC = [0 0 wBEC];

omega_CD = [0 0 wCD];

eqn11 = cross(omega_AB,B-A) + cross(omega_BEC,C-B) + cross(omega_CD,D-C) == 0;

loopSolution = solve(eqn11,[wBEC wCD]);

angularVelocity_BEC = double(loopSolution.wBEC)

angularVelocity_CD = double(loopSolution.wCD)

omegaBEC = [0 0 angularVelocity_BEC];

omegaCD = [0 0 angularVelocity_CD];


% Loop 2 DCEF​​G

syms wEF wFG

omega_EF = [0 0 wEF];

omega_FG = [0 0 wFG];

eqn12 = cross(omegaCD,C-D) + cross(omegaBEC,E-C) + cross(omega_EF,F-E) + cross(omega_FG,G-F) == 0;

loop2Solution = solve(eqn12,[wEF wFG]);

angularVelocity_EF = double(loop2Solution.wEF)

angularVelocity_FG = double(loop2Solution.wFG)

angVEL_EF = [0 0 angularVelocity_EF];

angVEL_FG = [0 0 angularVelocity_FG];


%% Angular Acceleration

% Loop 1 ABCDA

syms aBEC aCD

alpha_AB = [0 0 0];

alpha_BEC = [0 0 aBEC];

alpha_CD = [0 0 aCD];

a_B_A = cross(alpha_AB,B-A) + cross(omega_AB,cross(omega_AB,B-A));

a_C_B = cross(alpha_BEC,C-B) + cross(omegaBEC,cross(omegaBEC,C-B));

a_D_C = cross(alpha_CD,D-C) + cross(omegaCD,cross(omegaCD,D-C));

eqn13 = a_B_A + a_C_B + a_D_C == 0;

loop1AccSolution = solve(eqn13,[aBEC aCD]);

alphaBEC = double(loop1AccSolution.aBEC)

alphaCD = double(loop1AccSolution.aCD)

alphaBEC_vector = [0 0 alphaBEC];

alphaCD_vector = [0 0 alphaCD];


% Loop 2

syms aEF aFG

alpha_EF = [0 0 aEF];

alpha_FG = [0 0 aFG];

a_C_D = cross(alphaCD_vector,C-D) + cross(omegaCD,cross(omegaCD,C-D));

a_E_C = cross(alphaBEC_vector,E-C) + cross(omegaBEC,cross(omegaBEC,E-C));

a_F_E = cross(alpha_EF,F-E) + cross(angVEL_EF,cross(angVEL_EF,F-E));

a_G_F = cross(alpha_FG,G-F) + cross(angVEL_FG,cross(angVEL_FG,G-F));

eqn14 = a_C_D + a_E_C + a_F_E + a_G_F == 0;

loop2AccSolution = solve(eqn14,[aEF aFG]);

alphaEF = double(loop2AccSolution.aEF)

alphaFG = double(loop2AccSolution.aFG)

alphaEF_vector = [0 0 alphaEF];

alphaFG_vector = [0 0 alphaFG];


%% Velocity at Joint


% Velocity at Joint B

vB_A = cross(omega_AB,B-A)


% Velocity at Joint C

vC_D = cross(omegaCD,C-D)


% Velocity at Joint E

% vE_A = V_E_B + V_B_A

v_E_B = cross(omegaBEC,E-B);

vE_A = v_E_B + vB_A


% Velocity at Joint F

V_F_G = cross(angVEL_FG,F-G)


%% Velocity at Center of Mass


% Velocity of S1

vB_A = cross(omega_AB,B-A);

% vS1_A = V_S1_B + V_B_A

V_S1_B = cross(omega_AB,S1-B);

vS1_A = V_S1_B + vB_A


% Velocity of S2

vB_A = cross(omega_AB,B-A);

% vS2_A = V_S2_B + V_B_A

V_S2_B = cross(omegaBEC,S2-B);

vS2_A = V_S2_B + vB_A


% Velocity of S3

vC_D = cross(omegaCD,C-D);

% vS3_D = V_S3_C + V_C_D

V_S3_C = cross(omegaCD,S3-C);

vS3_D = V_S3_C + vC_D


% Velocity of S4

V_F_G = cross(angVEL_FG,F-G);

% vS4_G = V_S4_F + V_F_G

V_S4_F = cross(angVEL_EF,S4-F);

vS4_G = V_S4_F + V_F_G


% Velocity of S5

V_F_G = cross(angVEL_FG,F-G);

% vS5_G = V_S5_F + V_F_G

V_S5_F = cross(angVEL_FG,S5-F);

vS5_G = V_S5_F + V_F_G


%% Acceleration at Joint


% Acceleration at Joint B

aB_A = cross(alpha_AB,B-A) + cross(omega_AB,cross(omega_AB,B-A))


% Acceleration at Joint C

aC_D = cross(alphaCD_vector,C-D) + cross(omegaCD,cross(omegaCD,C-D))


% Acceleration at Joint E

% aE_A = a_E_B + a_B_A

a_E_B = cross(alphaBEC_vector,E-B) + cross(omegaBEC,cross(omegaBEC,E-B));

aE_A = a_E_B + aB_A


% Acceleration at Joint F

aF_G = cross(alphaFG_vector,F-G) + cross(angVEL_FG,cross(angVEL_FG,F-G))


%% Acceleration at Center of Mass


% Acceleration of S1

aB_A = cross(alpha_AB,B-A) + cross(omega_AB,cross(omega_AB,B-A));

% aS1_A = a_S1_B + a_B_A

a_S1_B = cross(alpha_AB,S1-B) + cross(omega_AB,cross(omega_AB,S1-B));

aS1_A = a_S1_B + aB_A


% Acceleration of S2

aB_A = cross(alpha_AB,B-A) + cross(omega_AB,cross(omega_AB,B-A));

% aS2_A = a_S2_B + a_B_A

a_S2_B = cross(alphaBEC_vector,S2-B) + cross(omegaBEC,cross(omegaBEC,S2-B));

aS2_A = a_S2_B + aB_A


% Acceleration of S3

aC_D = cross(alphaCD_vector,C-D) + cross(omegaCD,cross(omegaCD,C-D));

% aS3_D = a_S3_C + a_C_D

a_S3_C = cross(alphaCD_vector,S3-C) + cross(omegaCD,cross(omegaCD,S3-C));

aS3_D = a_S3_C + aC_D


% Acceleration of S4

aF_G = cross(alphaFG_vector,F-G) + cross(angVEL_FG,cross(angVEL_FG,F-G));

% aS4_G = a_S4_F + a_F_G

a_S4_F = cross(alphaEF_vector,S4-F) + cross(angVEL_EF,cross(angVEL_EF,S4-F));

aS4_G = a_S4_F + aF_G


% Acceleration of S5

aF_G = cross(alphaFG_vector,F-G) + cross(angVEL_FG,cross(angVEL_FG,F-G));

% aS5_G = a_S5_F + a_F_G

a_S5_F = cross(alphaFG_vector,S5-F) + cross(angVEL_FG,cross(angVEL_FG,S5-F));

aS5_G = a_S5_F + aF_G


%% Newton's Second Law / Dynamic Force Analysis

% Mass of Each Link

MassAB = 1;
MassBEC = 1;
MassCD = 1;
MassEF = 1;
MassFG = 1;

% Mass Moment of Inertia of Each Link

J_AB = 1;
J_BEC = 1;
J_CD = 1;
J_EF = 1;
J_FG = 1;

% Define Unknown Dynamic Forces and Input Torque

syms NFAx NFAy NFBx NFBy NFCx NFCy NFDx NFDy ...
     NFEx NFEy NFFx NFFy NFGx NFGy NTin

NForceA = [NFAx NFAy 0];
NForceB = [NFBx NFBy 0];
NForceC = [NFCx NFCy 0];
NForceD = [NFDx NFDy 0];
NForceE = [NFEx NFEy 0];
NForceF = [NFFx NFFy 0];
NForceG = [NFGx NFGy 0];

NInputTorque = [0 0 NTin];


%% Newton's Second Law Equations

% Link AB
% Sum of Forces = Mass * Acceleration of Center of Mass

eqn15 = NForceA + NForceB + WAB == MassAB * aS1_A;

% Sum of Moments about Center of Mass = J * alpha

eqn16 = cross(A-S1,NForceA) + ...
        cross(B-S1,NForceB) + ...
        NInputTorque == J_AB * alpha_AB;


% Link BEC

eqn17 = -NForceB + NForceC + NForceE + WBEC == ...
        MassBEC * aS2_A;

eqn18 = cross(B-S2,-NForceB) + ...
        cross(C-S2,NForceC) + ...
        cross(E-S2,NForceE) == ...
        J_BEC * alphaBEC_vector;


% Link CD

eqn19 = -NForceC + NForceD + WCD == ...
        MassCD * aS3_D;

eqn20 = cross(C-S3,-NForceC) + ...
        cross(D-S3,NForceD) == ...
        J_CD * alphaCD_vector;


% Link EF

eqn21 = -NForceE + NForceF + WEF == ...
        MassEF * aS4_G;

eqn22 = cross(E-S4,-NForceE) + ...
        cross(F-S4,NForceF) == ...
        J_EF * alphaEF_vector;


% Link FG

eqn23 = -NForceF + NForceG + WFG + AppliedForce == ...
        MassFG * aS5_G;

eqn24 = cross(F-S5,-NForceF) + ...
        cross(G-S5,NForceG) == ...
        J_FG * alphaFG_vector;


%% Collect Equations and Solve

NeqnMatrix = [eqn15, eqn16, ...
              eqn17, eqn18, ...
              eqn19, eqn20, ...
              eqn21, eqn22, ...
              eqn23, eqn24];

DynamicSolution = solve(NeqnMatrix, ...
    [NFAx NFAy NFBx NFBy NFCx NFCy NFDx NFDy ...
     NFEx NFEy NFFx NFFy NFGx NFGy NTin]);


%% Convert Solutions to Numeric Values

NForce_Ax = double(DynamicSolution.NFAx)
NForce_Ay = double(DynamicSolution.NFAy)

NForce_Bx = double(DynamicSolution.NFBx)
NForce_By = double(DynamicSolution.NFBy)

NForce_Cx = double(DynamicSolution.NFCx)
NForce_Cy = double(DynamicSolution.NFCy)

NForce_Dx = double(DynamicSolution.NFDx)
NForce_Dy = double(DynamicSolution.NFDy)

NForce_Ex = double(DynamicSolution.NFEx)
NForce_Ey = double(DynamicSolution.NFEy)

NForce_Fx = double(DynamicSolution.NFFx)
NForce_Fy = double(DynamicSolution.NFFy)

NForce_Gx = double(DynamicSolution.NFGx)
NForce_Gy = double(DynamicSolution.NFGy)

NInput_Torque = double(DynamicSolution.NTin)

initial_theta = atan2(B(2)-A(2),B(1)-A(1));

if (initial_theta<0)
    inputAngle = 2*pi + initial_theta;
else

    inputAngle = initial_theta;
end


for theta=1:1:360
    % new pos of joint B
    B_new = A + [lAB*cos(inputAngle+deg2rad(theta)) lAB*sin(inputAngle+deg2rad(theta)) 0]
    % new pos of joint C
    [Cx, Cy] = circcirc(B_new(1),B_new(2),lBC,D(1),D(2),lCD);

%% 
% checking if there is a NaN
circIntersect_x = any(isnan(vpa(Cx)));
circIntersect_y = any(isnan(vpa(Cy)));

if circIntersect_x==0 && circIntersect_y==0
    C_1 = [Cx(1) Cy(1) 0];
    C_2 = [Cy(2) Cy(2) 0];

    % distance to determine whether C_1 or C_2 is correct
    dist1 = norm(C_1-C);
    dist2 = norm(C_2-C);

    if(dist1<dist2)
        C_new = vpa(C_1);
    else
        C_new = vpa(C_2);
    end

    % New position of Joint E using B_new and C_new

    [Ex,Ey] = circcirc(B_new(1),B_new(2),lBE,C_new(1),C_new(2),lCE);

    % checking if there is a NaN
    circIntersect_x = any(isnan(vpa(Ex)));
    circIntersect_y = any(isnan(vpa(Ey)));

    if circIntersect_x==0 && circIntersect_y==0
        E_1 = [Ex(1) Ey(1) 0];
        E_2 = [Ex(2) Ey(2) 0];

        % distance to determine whether C_1 or C_2 is correct
        dist1 = norm(E_1-E);
        dist2 = norm(E_2-E);

        if(dist1<dist2)
            E_new = vpa(E_1);
        else
            E_new = vpa(E_2);
        end

    % New position of Joint F using E_new and G_new

    [Fx,Fy] = circcirc(E_new(1),E_new(2),lEF,G(1),G(2),lFG);

    % checking if there is a NaN
    circIntersect_x = any(isnan(vpa(Fx)));
    circIntersect_y = any(isnan(vpa(Fy)));

    if circIntersect_x==0 && circIntersect_y==0
        F_1 = [Fx(1) Fy(1) 0];
        F_2 = [Fx(2) Fy(2) 0];

        % distance to determine whether C_1 or C_2 is correct
        dist1 = norm(F_1-F);
        dist2 = norm(F_2-F);

        if(dist1<dist2)
            F_new = vpa(F_1);
        else
            F_new = vpa(F_2);
        end

        %Store value for plotting

        new_B_x(theta+1) = B_new(1);
        new_B_y(theta+1) = B_new(2);
        new_C_x(theta+1) = C_new(1);
        new_C_y(theta+1) = C_new(2);
        new_E_x(theta+1) = E_new(1);
        new_E_y(theta+1) = E_new(2);
        new_F_x(theta+1) = F_new(1);
        new_F_y(theta+1) = F_new(2);

        B=B_new;
        C=C_new;
        E=E_new;
        F=F_new;

        %Static Equilibrium Code

        eqn1 = ForceA + ForceB + WAB == 0;

        eqn2 = cross(A-S1,ForceA) + cross(B-S1,ForceB) + InputTorque == 0;

        eqn3 = -ForceB + ForceC + ForceE + WBEC == 0;

        eqn4 = cross(B-S2,-ForceB) + cross(C-S2,ForceC) + cross(E-S2,ForceE) == 0;

        eqn5 = -ForceC + ForceD + WCD == 0;

        eqn6 = cross(C-S3,-ForceC) + cross(D-S3,ForceD) == 0;

        eqn7 = -ForceE + ForceF + WEF == 0;

        eqn8 = cross(E-S4,-ForceE) + cross(F-S4,ForceF) == 0;

        eqn9 = -ForceF + ForceG + WFG + AppliedForce == 0;

        eqn10 = cross(F-S5,-ForceF) + cross(G-S5,ForceG) == 0;

        eqnMatrix = [eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8, eqn9, eqn10];

        StaticSolution = solve(eqnMatrix, [FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin]);

        Force_Ax(theta+1) = double(StaticSolution.FAx);
        Force_Ay(theta+1) = double(StaticSolution.FAy);

        Force_Bx(theta+1) = double(StaticSolution.FBx);
        Force_By(theta+1) = double(StaticSolution.FBy);

        Force_Cx(theta+1) = double(StaticSolution.FCx);
        Force_Cy(theta+1) = double(StaticSolution.FCy);

        Force_Dx(theta+1) = double(StaticSolution.FDx);
        Force_Dy(theta+1) = double(StaticSolution.FDy);

        Force_Ex(theta+1) = double(StaticSolution.FEx);
        Force_Ey(theta+1) = double(StaticSolution.FEy);

        Force_Fx(theta+1) = double(StaticSolution.FFx);
        Force_Fy(theta+1) = double(StaticSolution.FFy);

        Force_Gx(theta+1) = double(StaticSolution.FGx);
        Force_Gy(theta+1) = double(StaticSolution.FGy);

        Input_Torque = double(StaticSolution.Tin);

        %Velocity and Acc Equation

        %Newton's Second Law


    else

        max_theta = theta - 1;

        fprintf('\nF cannot be determined at theta = %d degrees\n',theta);
        fprintf('Maximum valid input rotation = %d degrees\n',max_theta);

        break
    end


else

    max_theta = theta - 1;

    fprintf('\nE cannot be determined at theta = %d degrees\n',theta);
    fprintf('Maximum valid input rotation = %d degrees\n',max_theta);

    break
end


else

    max_theta = theta - 1;

    fprintf('\nC cannot be determined at theta = %d degrees\n',theta);
    fprintf('Maximum valid input rotation = %d degrees\n',max_theta);

    break
end


end




