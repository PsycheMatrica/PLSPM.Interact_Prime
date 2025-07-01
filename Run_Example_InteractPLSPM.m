%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Illustration for PLSPM.Interact_Prime package                           %
%   Author: Gyeongcheol Cho & Heungsun Hwang                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:                                                            %
%  - This code aims to illustrate how to use PLSPM.Interact_Prime package.%
%  - The dataset is generated from the model used in Shen, Cho, and Hwang %
%     (2025).                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References                                                              %
%     * Chin, W. W., Marcolin, B. L., & Newsted, P. R. (2003). A Partial  %
%         Least Squares Latent Variable Modeling Approach for Measuring   %
%         Interaction Effects: Results from a Monte Carlo Simulation Study% 
%         and an Electronic-Mail Emotion/Adoption Study. Information      %
%         Systems Research, 14(2), 189–217.                               %
%         https://doi.org/10.1287/isre.14.2.189.16018                     %
%     * Shen, Z., Cho, G., & Hwang, H. (2025). Comparison of              %
%       Component-Based Structural Equation Modeling Methods in Testing   %
%       Component Interaction Effects. Structural Equation Modeling: A    %
%       Multidisciplinary Journal, 1–13.                                  %
%       https://doi.org/10.1080/10705511.2025.2497088                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
z0=readmatrix('Example_data.csv');
P = 3;
P_int = 1;
P_total = 4;
Jp = 5;
J = Jp * P;

wp  = ones(Jp,1);
W0 = blkdiag(wp,wp,wp)*99;
C0 = W0';
nnlv_index=[1,2];
B0 = zeros(P_total,P);
B0([1,2,4],3)=99;
ind_sign=[1,6,11];
N_Boot = 100;
modetype=ones(1,P);
scheme=3;
Flag_Parallel=false;

Max_iter = 100; 
Min_limit = .00001;

Results = InteractPLSPM(z0, W0, B0, modetype, scheme, nnlv_index, ind_sign, N_Boot, Max_iter, Min_limit, Flag_Parallel);
INI=Results.INI;
TABLE=Results.TABLE;
ETC=Results.ETC;
 
INI.W 
INI.C
INI.B 

TABLE.W
TABLE.C
TABLE.B