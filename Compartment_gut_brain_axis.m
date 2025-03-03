function dy = Compartment_gut_brain_axis(time, y, sympathetic)

% Initialize variables for Pyloric Sphincter (Phasic and Tonic) 
% These are not state variables, but auxiliary variables that may be used in calculations.
% Since antrum phasic contractions does not depend on any factors, we set the baseline stimilus as 0
f_i_PS_phasic = 0; 
f_e_PS_phasic = 0;

% for parasympethetic
f_e_PS_tonic = 0;

% -------------------- Fundus (State Variables) --------------------
% These are state variables representing different biochemical interactions
% related to calcium-Calmodulin (CaM) and myosin regulation in the fundus.
Ca2_CaM_comp_1 = y(1, :);   % Ca2+-bound CaM
CaM_comp_1 = y(2, :);       % Free CaM
Ca4_CaM_comp_1 = y(3, :);   % Ca4+-bound CaM
CaM_MLCK_comp_1 = y(4, :);  % CaM bound to MLCK
MLCK_comp_1 = y(5, :);      % Free MLCK
Ca2_CaM_MLCK_comp_1 = y(6, :);  % Ca2+-CaM-MLCK complex
Ca4_CaM_MLCK_comp_1 = y(7, :);  % Ca4+-CaM-MLCK complex
CaM_Buf_comp_1 = y(8, :);   % CaM bound to buffer
Buf_comp_1 = y(9, :);       % Free buffer
Mp_comp_1 = y(10, :);       % Phosphorylated myosin
M_comp_1 = y(11, :);        % Unphosphorylated myosin
AMp_comp_1 = y(12, :);      % Actomyosin complex (phosphorylated)
AM_comp_1 = y(13, :);       % Actomyosin complex (unphosphorylated)

% -------------------- Antrum (State Variables) --------------------
% These state variables define membrane potentials and calcium signaling in 
% the interstitial cells of Cajal (ICC) and smooth muscle cells (SMC) of the antrum.
V_m_ICC_comp_2 = y(14, :);  % ICC membrane potential
V_m_SMC_comp_2 = y(15, :);  % SMC membrane potential
Ca_i_SM_comp_2 = y(16, :);  % Intracellular Ca2+ in smooth muscle
f_Ltype_SM_comp_2 = y(17, :); % L-type Ca2+ channel (fast activation)
d_Ltype_SM_comp_2 = y(18, :); % L-type Ca2+ channel (slow activation)
f_ca_Ltype_SM_comp_2 = y(19, :); % L-type Ca2+ channel (Ca-dependent)
f_LVA_SM_comp_2 = y(20, :);  % Low-voltage-activated Ca2+ channel (fast)
d_LVA_SM_comp_2 = y(21, :);  % Low-voltage-activated Ca2+ channel (slow)

% State variables for calcium-Calmodulin interactions in the antrum
Ca2_CaM_comp_2 = y(22, :);
CaM_comp_2 = y(23, :);
Ca4_CaM_comp_2 = y(24, :);
CaM_MLCK_comp_2 = y(25, :);
MLCK_comp_2 = y(26, :);
Ca2_CaM_MLCK_comp_2 = y(27, :);
Ca4_CaM_MLCK_comp_2 = y(28, :);
CaM_Buf_comp_2 = y(29, :);
Buf_comp_2 = y(30, :);
Mp_comp_2 = y(31, :);
M_comp_2 = y(32, :);
AMp_comp_2 = y(33, :);
AM_comp_2 = y(34, :);
stretch_comp_2 = y(35, :);  % Stretch component in antrum

% -------------------- Pyloric Sphincter (Tonic) (State Variables) --------------------
% State variables representing CaM interactions in the tonic region of the pyloric sphincter.
Ca2_CaM_tonic_comp_3 = y(36, :);
CaM_tonic_comp_3 = y(37, :);
Ca4_CaM_tonic_comp_3 = y(38, :);
CaM_MLCK_tonic_comp_3 = y(39, :);
MLCK_tonic_comp_3 = y(40, :);
Ca2_CaM_MLCK_tonic_comp_3 = y(41, :);
Ca4_CaM_MLCK_tonic_comp_3 = y(42, :);
CaM_Buf_tonic_comp_3 = y(43, :);
Buf_tonic_comp_3 = y(44, :);
Mp_tonic_comp_3 = y(45, :);
M_tonic_comp_3 = y(46, :);
AMp_tonic_comp_3 = y(47, :);
AM_tonic_comp_3 = y(48, :);

% -------------------- Pyloric Sphincter (Phasic) (State Variables) --------------------
% State variables representing electrical and biochemical dynamics in the phasic region.
V_m_ICC_phasic_comp_3 = y(49, :); 
V_m_SMC_phasic_comp_3 = y(50, :);
Ca_i_SM_phasic_comp_3 = y(51, :);
f_Ltype_SM_phasic_comp_3 = y(52, :);
d_Ltype_SM_phasic_comp_3 = y(53, :);
f_ca_Ltype_SM_phasic_comp_3 = y(54, :);
f_LVA_SM_phasic_comp_3 = y(55, :);
d_LVA_SM_phasic_comp_3 = y(56, :);

% State variables for calcium-Calmodulin interactions in the phasic region of the pyloric sphincter
Ca2_CaM_phasic_comp_3 = y(57, :);
CaM_phasic_comp_3 = y(58, :);
Ca4_CaM_phasic_comp_3 = y(59, :);
CaM_MLCK_phasic_comp_3 = y(60, :);
MLCK_phasic_comp_3 = y(61, :);
Ca2_CaM_MLCK_phasic_comp_3 = y(62, :);
Ca4_CaM_MLCK_phasic_comp_3 = y(63, :);
CaM_Buf_phasic_comp_3 = y(64, :);
Buf_phasic_comp_3 = y(65, :);
Mp_phasic_comp_3 = y(66, :);
M_phasic_comp_3 = y(67, :);
AMp_phasic_comp_3 = y(68, :);
AM_phasic_comp_3 = y(69, :);
stretch_phasic_comp_3 = y(70, :); % Stretch component in phasic region

Volume_fluid = y(71, :); % Volume of fluid component

% -------------------- Volume Calculation --------------------
% Total volume calculation
% The previously commented-out equation suggests that Volume_total could depend on Volume_fluid.
% However, it is now assigned a fixed value of 600, overriding any dependency on Volume_fluid.
%Volume_total = 600;


%Dynamic value
V_gas = 210; %mL

Volume_total = Volume_fluid + V_gas;


%% Constants

% Calcium concentration SMC
F=96484.6;%Funits(inEnvironment)
Vol_SM=3.5e-3;%volumeunits(inMembrane)
tau_f_Ltype_SM=225.6231e-3;%timeunits(infLtypeSM)
tau_d_Ltype_SM=1.2331e-3;%timeunits(indLtypeSM)
tau_f_ca_Ltype_SM=5.247e-3;%timeunits(infcaLtypeSM)
tau_d_LVA_SM=7.8706e-3;%timeunits(indLVASM)
G_max_Ltype_2=65.0;%conductanceunits(GmaxLtypeinILtypeSM)
RToF=26.7137;%voltageunits(inEnvironment)
T_correction_Ca=2.6235;%dimensionless(inEnvironment)
Ca_o_SM=2.5;%millimolar(inEnvironment)
G_max_LVA=0.18;%conductanceunits(inILVASM)
J_max_CaSR=317.05;%millimolarpersecond(inJCaSRSM)

%baseline frequency for rhythmic activity in stomach
base_cpm = 3;

%baseline amplitude for rhythmic activity in stomach
base_amplitude_ICC = 48;

%% Mechanical activity parameters

% --- Reactions (MLC activation) Table 1 of Gajendiran & Buist --- //

 % CaM_comp_1 + 2Ca <-> Ca2-CaM_comp_1
  Kf_Ca2CaM = 12.0;                    % uM/s
  Kr_Ca2CaM = 12.0;                    % 1/s
  Kf_Ca2Ca4 = 480.0;                   % uM/s
  Kr_Ca2Ca4 = 1200.0;                  % 1/s
  Kf_CaMMLCK = 5.0;                    % uM/s
  Kr_CaMMLCK = 135.0;                  % 1/s
  Kf_Ca2CaMtoCa2_CaM_MLCK = 840.0;     % uM/s
  Kr_Ca2CaMtoCa2_CaM_MLCK = 45.4;      % 1/s
  Kf_MLCKact = 28.0;                   % uM/s
  Kr_MLCKact = 0.0308;                 % 1/s
  Kf_Ca2CaMMLCK = 120.0;               % uM/s
  Kr_Ca2CaMMLCK = 4.0;                 % 1/s
  Kf_CaCaMMLCK = 7.5;                  % uM/s
  Kr_CaCaMMLCK = 3.75;                 % 1/s
  Kf_CaMBuf = 5.0;                     % uM/s
  Kr_CaMBuf = 25.0;                    % 1/s
  Kf_CaMBufdis = 7.6;                  % uM/s
  Kr_CaMBufdis = 22.8;                 % 1/s
  Vmax_MLCKact = 27.0;                 % 1/s
  Km_MLCKact = 10.0;                   % uM
  HM_k3 = 15.0;                        % 1/s
  HM_k4 = 5.0;                         % 1/s
  Vmax_MLCP = 16.0;                    % 1/s
  Km_MLCP = 15.0;                      % uM
  HM_k7 = 10.0;                        % 1/s
  HM_k8 = 0.0;                         % 1/s
%% Afferent - Fundus

% Enforce a minimum gastric volume of 80
if Volume_total <= 80
    Volume_total = 80;
end

% -------------------- NANC Volume Coefficients --------------------
% Coefficients for polynomial approximations of fundus function (NANC volume)
% These define different equations for low and high gastric volumes.

% Coefficients for high gastric volume condition
k_vol_high_1 = 1.98824066882335E-17;
k_vol_high_2 = -1.94240191873764E-13;
k_vol_high_3 = 8.43067475381809E-10;
k_vol_high_4 = -2.13371855645172E-06;
k_vol_high_5 = 0.00347024704364;
k_vol_high_6 = -3.76119327725484;
k_vol_high_7 = 2716.64090413506;
k_vol_high_8 = -1260913.39204722;
k_vol_high_9 = 341260848.375847;
k_vol_high_10 = -41033260585.2664;

% Coefficients for low gastric volume condition
k_vol_low_1 = -5.27526938512081E-22;
k_vol_low_2 = 2.72372868184093E-18;
k_vol_low_3 = -5.88796879264595E-15;
k_vol_low_4 = 6.94303089718836E-12;
k_vol_low_5 = -4.83821602296794E-09;
k_vol_low_6 = 0.0000020068929203711;
k_vol_low_7 = -0.00046754318758084;
k_vol_low_8 = 0.058051209557502;
k_vol_low_9 = -2.48610055439772;

% -------------------- Fundus Function Calculation --------------------
% Computes f_i_fundus based on Volume_total using different polynomial models
% Low-volume equation used if Volume_total ≤ 997
% High-volume equation used otherwise
% Piecewise plunimial for NANC

if Volume_total <= 997
    f_i_fundus = k_vol_low_1 .* Volume_total.^8 + ...
                 k_vol_low_2 .* Volume_total.^7 + ...
                 k_vol_low_3 .* Volume_total.^6 + ...
                 k_vol_low_4 .* Volume_total.^5 + ...
                 k_vol_low_5 .* Volume_total.^4 + ...
                 k_vol_low_6 .* Volume_total.^3 + ...
                 k_vol_low_7 .* Volume_total.^2 + ...
                 k_vol_low_8 .* Volume_total + ...
                 k_vol_low_9;
else
    f_i_fundus = k_vol_high_1 .* Volume_total.^9 + ...
                 k_vol_high_2 .* Volume_total.^8 + ...
                 k_vol_high_3 .* Volume_total.^7 + ...
                 k_vol_high_4 .* Volume_total.^6 + ...
                 k_vol_high_5 .* Volume_total.^5 + ...
                 k_vol_high_6 .* Volume_total.^4 + ...
                 k_vol_high_7 .* Volume_total.^3 + ...
                 k_vol_high_8 .* Volume_total.^2 + ...
                 k_vol_high_9 .* Volume_total + ...
                 k_vol_high_10;
end

% Cholinergic fundus neuron stimulus (this is set to 0 as only the
% NANC pathway participates in fundus relaxation based on the assumption made in the paper)
f_e_fundus = 0;


%% Afferent - PS (Pyloric Sphincter)
% This section models the afferent chemical response in the pyloric sphincter.

% -------------------- Caloric Absorption --------------------
g_cal = 400./600; % Initial gastric caloric flow ratio
% g_cal = 5./600; % Alternative value (commented out)

gamma_max = 3; % Maximum gastric energy absorption rate (kcal/min)
omega = 1.2;   % Gastric emptying rate factor

% Ensure a minimum caloric flow to prevent division errors
if g_cal <= 0.0514299628120262
    g_cal = 0.0514299628120262;
end

% Compute caloric absorption rate
Q_cal = (gamma_max .* omega) ./ (g_cal .* 60); % Normalized caloric absorption rate

% Compute afferent chemical signaling based on caloric absorption
% Previous equation commented out
% f_afferent_chem = 30.626664109729 .* Q_cal - 5.73130556439374;
f_afferent_chem = 26.3941592 .* Q_cal - 0.7923526592644;

% Prevent negative afferent chemical response
if f_afferent_chem <= 0
    f_afferent_chem = 0;
end

% -------------------- Polynomial Approximation for Tonic Component --------------------
% Coefficients for polynomial function defining the tonic response
k_chem_flow_1 =  2.628606668248840e-10;
k_chem_flow_2 = -3.265010083604308e-08;
k_chem_flow_3 =  1.700068150165396e-06;
k_chem_flow_4 = -4.814655421754510e-05;
k_chem_flow_5 =  8.063641587893079e-04;
k_chem_flow_6 = -0.008150916917974;
k_chem_flow_7 =  0.049024001598177;
k_chem_flow_8 = -0.168738983217585;
k_chem_flow_9 =  0.398146785989247;
k_chem_flow_10 =  0.021430144593178;

% Compute tonic afferent response using a 9th-degree polynomial function
f_i_PS_tonic = k_chem_flow_1 .* f_afferent_chem.^9 + ...
               k_chem_flow_2 .* f_afferent_chem.^8 + ...
               k_chem_flow_3 .* f_afferent_chem.^7 + ...
               k_chem_flow_4 .* f_afferent_chem.^6 + ...
               k_chem_flow_5 .* f_afferent_chem.^5 + ...
               k_chem_flow_6 .* f_afferent_chem.^4 + ...
               k_chem_flow_7 .* f_afferent_chem.^3 + ...
               k_chem_flow_8 .* f_afferent_chem.^2 + ...
               k_chem_flow_9 .* f_afferent_chem + ...
               k_chem_flow_10;

% Enforce a minimum threshold for tonic response
if f_i_PS_tonic <= 0
    f_i_PS_tonic = 0;
end

% -------------------- Sympathetic Nervous System Influence --------------------
% When the sympathetic nervous system is active, afferent tonic input is suppressed.
if sympathetic == 1
    f_i_PS_tonic = 0;  % No afferent tonic response
    f_e_PS_tonic = 10; % Increased efferent output
end


%% Efferent Neuron Signaling - Fundus
% Models efferent signaling in the fundus, including Ach, VIP, NO, cGMP, and MLCP.

% -------------------- Acetylcholine (Ach) Release --------------------
% Efferent fundus signaling parameters for Ach release
a_fundus_ach_tonic = 1.03130949312859e3;
b_fundus_ach_tonic = 0.165129489206726e3;
n_fundus_ach_tonic = 0.00306813263509e3;

% Compute Ach release from efferent neurons
ach_fundus_tonic = (a_fundus_ach_tonic .* power(f_e_fundus, n_fundus_ach_tonic)) ./ ...
                    (power(f_e_fundus, n_fundus_ach_tonic) + power(b_fundus_ach_tonic, n_fundus_ach_tonic));

% -------------------- Calcium (Ca) Release Mediated by Ach --------------------
% Parameters for calcium release in response to Ach
a_fundus_Ca_tonic = 0.000233146317783e5;
b_fundus_Ca_tonic = 4.893974398e5;
n_fundus_Ca_tonic = 0.346091;

% Compute Ca release in response to Ach
Ca_fundus_tonic = (a_fundus_Ca_tonic .* power(ach_fundus_tonic, n_fundus_Ca_tonic)) ./ ...
                   (power(ach_fundus_tonic, n_fundus_Ca_tonic) + power(b_fundus_Ca_tonic, n_fundus_Ca_tonic));

% -------------------- VIP (Vasoactive Intestinal Peptide) Release --------------------
% VIP release parameters based on afferent input
a_VIP_fundus = 67.372030641666839;
b_VIP_fundus = 19.896399126795536;
n_VIP_fundus = 2.353830304145993;


% Compute VIP release in response to afferent signaling
VIP_fundus_tonic = (a_VIP_fundus .* power(f_i_fundus, n_VIP_fundus)) ./ ...
                    (power(b_VIP_fundus, n_VIP_fundus) + power(f_i_fundus, n_VIP_fundus));



% -------------------- VIP Facilitation and Inhibition --------------------
% Parameters for VIP facilitation
n_VIP_fac_fundus = 0.824092403711139;
b_VIP_fac_fundus = 1.254001767733051;


% Compute the VIP-mediated facilitation factor
factor_VIP_fundus = power(VIP_fundus_tonic, n_VIP_fac_fundus) ./ ...
                    (power(VIP_fundus_tonic, n_VIP_fac_fundus) + power(b_VIP_fac_fundus, n_VIP_fac_fundus));

% Compute the inhibitory effect of VIP
factor_inhibitory = 1 - factor_VIP_fundus;

% -------------------- Nitric Oxide (NO) Computation --------------------
% Parameters for NO production in the fundus (in nM)
a_NO_fundus = 0.667664813;
b_NO_fundus = 18;
n_NO_fundus = 1.501113355;


% Compute NO concentration based on afferent input
NO_fundus_tonic_nM = (a_NO_fundus .* power(f_i_fundus, n_NO_fundus)) ./ ...
                      (power(b_NO_fundus, n_NO_fundus) + power(f_i_fundus, n_NO_fundus));

% Convert NO from uM to nM
NO_fundus_tonic = NO_fundus_tonic_nM * 1000;

% -------------------- cGMP Production from NO --------------------
% Parameters for cGMP production in response to NO
a_cGMP_fundus = 26.5732;
b_cGMP_fundus = 21.6111;
n_cGMP_fundus = 1.0022;

% Compute cGMP production in response to NO
cGMP_fundus_tonic = (a_cGMP_fundus * power(NO_fundus_tonic, n_cGMP_fundus)) ./ ...
                     (power(NO_fundus_tonic, n_cGMP_fundus) + power(b_cGMP_fundus, n_cGMP_fundus));

% -------------------- MLCP Activation --------------------
% Parameters for MLCP activation
n_mlcp = 2;
K_m_mlcp = 5.5;


% Compute MLCP activation fraction based on cGMP
R_mlcp_fundus_tonic = (power(cGMP_fundus_tonic, n_mlcp)) ./ ...
                        (power(cGMP_fundus_tonic, n_mlcp) + power(K_m_mlcp, n_mlcp));  

% Compute MLCP concentration in the fundus
MLCP_comp_1 = 7.5 + 7.5 * R_mlcp_fundus_tonic;




%% Efferent Neuron Signaling - Pyloric Sphincter
% This section models the effects of Acetylcholine (Ach), NANC tonic inhibitory, NO, and Purinergic in the pyloric sphincter.

% -------------------- Tonic Component --------------------

% Parameters for Ach release in the tonic region of the pyloric sphincter
a_f_e_PS_tonic = 601.1770659194432;
b_f_e_PS_tonic = 379.0511541958966;
n_f_e_PS_tonic = 2.1129902525713;

% Ensure f_e_PS_tonic is nonzero to prevent division errors
f_e_PS_tonic = max(f_e_PS_tonic, 1e-6);

% Compute tonic Ach release based on efferent input
Ach_PS_tonic = (a_f_e_PS_tonic .* power(f_e_PS_tonic, n_f_e_PS_tonic)) ./ ...
               (power(f_e_PS_tonic, n_f_e_PS_tonic) + power(b_f_e_PS_tonic, n_f_e_PS_tonic));

% Parameters for NANC modulation
a_NANC_PS_tonic = 1;
b_NANC_PS_tonic = 1.027745153459091;
n_NANC_PS_tonic = 1.569900559816452;

% Ensure f_i_PS_tonic is nonzero to prevent division errors
f_i_PS_tonic = max(f_i_PS_tonic, 1e-6);

% Compute NANC factor for tonic input
factor_NANC_PS_tonic = (a_NANC_PS_tonic .* power(f_i_PS_tonic, n_NANC_PS_tonic)) ./ ...
                        (power(f_i_PS_tonic, n_NANC_PS_tonic) + power(b_NANC_PS_tonic, n_NANC_PS_tonic));

% -------------------- Phasic Component --------------------

% Parameters for NO production in the phasic region
a_NO_PS = 0.667664813;
b_NO_PS = 43.82034347;
n_NO_PS = 1.501113355;


% Compute NO release in phasic region
NO_PS_phasic_uM = (a_NO_PS .* power(f_i_PS_phasic, n_NO_PS)) ./ ...
                   (power(b_NO_PS, n_NO_PS) + power(f_i_PS_phasic, n_NO_PS));

% Convert NO from uM to nM
NO_PS_phasic = NO_PS_phasic_uM * 1000;


% Parameters for cGMP production in response to NO
a_cGMP_PS = 26.5732;
b_cGMP_PS = 21.6111;
n_cGMP_PS = 1.0022;


% Compute cGMP in response to NO
cGMP_PS = (a_cGMP_PS * power(NO_PS_phasic, n_cGMP_PS)) ./ ...
           (power(NO_PS_phasic, n_cGMP_PS) + power(b_cGMP_PS, n_cGMP_PS));

% Compute MLCP activation in the phasic region based on cGMP levels
R_mlcp_PS_phasic = (power(cGMP_PS, n_mlcp)) / ...
                   (power(cGMP_PS, n_mlcp) + power(K_m_mlcp, n_mlcp));  

% Compute MLCP concentration in the phasic component of the pyloric sphincter
MLCP_phasic_comp_3 = 7.5 + 2.5 * R_mlcp_PS_phasic;

% Parameters for NO production in the phasic region (amplitude and frequency ICC)

a_i_amplitude_ICC_PS = 0.385470815230167;
b_i_amplitude_ICC_PS = 0.037915896;
n_i_amplitude_ICC_PS = 0.74976983;

a_i_cpm_PS = 0.1228;
b_i_cpm_PS = 3.5258;
n_i_cpm_PS = 1.5264;

% Compute inhibitory modulation of ICC amplitude by NO in the phasic region
phi_i_amplitude_ICC_PS = (a_i_amplitude_ICC_PS .* power(NO_PS_phasic_uM, n_i_amplitude_ICC_PS)) ./ ...
                          (power(b_i_amplitude_ICC_PS, n_i_amplitude_ICC_PS) + power(NO_PS_phasic_uM, n_i_amplitude_ICC_PS));

% Compute inhibitory modulation of ICC frequency by NO in the phasic region
phi_i_cpm_PS = (a_i_cpm_PS .* power(NO_PS_phasic_uM, n_i_cpm_PS)) ./ ...
               (power(b_i_cpm_PS, n_i_cpm_PS) + power(NO_PS_phasic_uM, n_i_cpm_PS));

% Parameters for purinergic signaling in the phasic component
a_Pur_PS = 0.525315108793863;
b_Pur_PS = 8.76385759972606;
n_Pur_PS = 1.74061375397748;

% Parameters for gap junction inhibition in smooth muscle cells (SMCs)
a_i_gap_SMC_PS = 167.5930670258801;
b_i_gap_SMC_PS = 907.9276538870989;
n_i_gap_SMC_PS = 0.8399709278706;

% Compute purinergic signaling influence on the phasic component
Pur_PS = (a_Pur_PS .* power(f_i_PS_phasic, n_Pur_PS)) ./ ...
         (power(b_Pur_PS, n_Pur_PS) + power(f_i_PS_phasic, n_Pur_PS));

% Compute gap junction inhibition effect on smooth muscle cells in the phasic component
phi_i_gap_SMC_PS = (a_i_gap_SMC_PS .* power(Pur_PS, n_i_gap_SMC_PS)) ./ ...
                    (power(b_i_gap_SMC_PS, n_i_gap_SMC_PS) + power(Pur_PS, n_i_gap_SMC_PS));

% Parameters for acetylcholine (Ach) modulation in the phasic region
a_ach_PS = 156.4317944;
b_ach_PS = 3.538564077;
n_ach_PS = 3.993133727;

% Parameters for ICC amplitude modulation by excitatory input
a_e_amplitude_ICC_PS = 0.15;
b_e_amplitude_ICC_PS = 0.574457529640854;
n_e_amplitude_ICC_PS = 0.33;

% Parameters for ICC frequency modulation by excitatory input
a_e_cpm_PS = 1.3;
b_e_cpm_PS = 96;
n_e_cpm_PS = 1;

% Compute Ach release in the phasic region based on efferent input
ach_PS_phasic = (a_ach_PS .* power(f_e_PS_phasic, n_ach_PS)) ./ ...
                 (power(b_ach_PS, n_ach_PS) + power(f_e_PS_phasic, n_ach_PS));

% Compute total Ach concentration, combining tonic and phasic components
ach_PS = ach_PS_phasic + Ach_PS_tonic;

% Parameters for Ach-induced calcium release in smooth muscle cells (SMCs)
a_ach_PS_tonic = 0.7091635609648;
b_ach_PS_tonic = 126.168746846009;
n_ach_PS_tonic = 0.2453486941887;

% Compute calcium release in SMCs in response to Ach
Ca_SMC_PS_ach_tonic = (a_ach_PS_tonic .* power(ach_PS, n_ach_PS_tonic)) ./ ...
                      (power(ach_PS, n_ach_PS_tonic) + power(b_ach_PS_tonic, n_ach_PS_tonic));

% Compute excitatory modulation of ICC amplitude by Ach
phi_e_amplitude_ICC_PS = (a_e_amplitude_ICC_PS .* power(ach_PS, n_e_amplitude_ICC_PS)) ./ ...
                          (power(b_e_amplitude_ICC_PS, n_e_amplitude_ICC_PS) + power(ach_PS, n_e_amplitude_ICC_PS));

% Compute excitatory modulation of ICC frequency by Ach
phi_e_cpm_PS = (a_e_cpm_PS .* power(ach_PS, n_e_cpm_PS)) ./ ...
               (power(b_e_cpm_PS, n_e_cpm_PS) + power(ach_PS, n_e_cpm_PS));

% Compute final ICC frequency (cycles per minute) incorporating excitatory and inhibitory effects
cpm_PS = base_cpm .* (phi_e_cpm_PS - phi_i_cpm_PS + 1);

% Compute final ICC amplitude incorporating excitatory and inhibitory effects
amplitude_ICC_PS = base_amplitude_ICC .* (phi_e_amplitude_ICC_PS - phi_i_amplitude_ICC_PS + 1);

% Compute inhibitory factor from gap junction modulation in SMCs
factor_i_gap_SMC_PS = (1 - phi_i_gap_SMC_PS);



%% Electrical Activity Parameters

% Define the cyclic time (in seconds) for different components
time_cyc_comp_2 = time;          % Cyclic time for component 2 (Antrum)
time_cyc_phasic_comp_3 = time;   % Cyclic time for phasic component 3 (Pyloric Sphincter)

% Time duration for channel opening (in seconds)
t_open = 7;  

% Decay constants for ICC and SMC electrical activity (unitless)
s_ICC = 0.05;  % Decay constant for interstitial cells of Cajal (ICC)
s_SMC = 0.05;  % Decay constant for smooth muscle cells (SMC)

% -------------------- Electrical Coupling (Gap Junctions) --------------------

% Electrical coupling conductance (unit: Siemens)
G_coup = 1.255;  % Assigning gap junction conductance

% -------------------- Resting Membrane Potentials --------------------

% Resting membrane potential of ICC and SMC (unit: mV)
V_m_rest_ICC = -70;  % Resting potential of ICC cells
V_m_rest_SMC = -70;  % Resting potential of SMC cells



%% Mechanical Activity Parameters

% -------------------- NLVM Model Parameters --------------------

% Parameters for component 2 (Antrum) and its interactions
S_3_comp_2_2_3_4 = 0.04;  % Scaling factor for parameter 3
S_4_comp_2_2_3_4 = 4.199;     % Scaling factor for parameter 4
S_5_comp_2_2_3_4 = -0.187;    % Scaling factor for parameter 5
S_6_comp_2_2_3_4 = 5;         % Scaling factor for parameter 6

L_t_comp_2 = 2.1;  % Tissue length in cm for component 2 (antrum)

% Parameters for phasic component 3 (Pyloric Sphincter, same as antrum) 
S_3_phasic_comp_3_2_3_4 = 0.04;
S_4_phasic_comp_3_2_3_4 = 4.199;
S_5_phasic_comp_3_2_3_4 = -0.187;
S_6_phasic_comp_3_2_3_4 = 5;

L_t_phasic_comp_3 = 6;  % Tissue length in cm for component 3 (PS)

% Duration of one phasic cycle based on contraction frequency
t_end_phasic_comp_3 = 60 ./ cpm_PS;  

% Start time for the phasic component cycle
t_start_phasic_comp_3 = 5;  

% -------------------- Phasic Contraction Cycle Logic --------------------

% If the cycle time is less than the start time, inhibit ICC activity
if (time_cyc_phasic_comp_3 < t_start_phasic_comp_3) 
    I_ICC_phasic_comp_3 = 0;
else 
    % Adjust the cycle time by subtracting the start time
    time_cyc_phasic_comp_3 = time_cyc_phasic_comp_3 - t_start_phasic_comp_3;

    % Restart the wave cycle after `t_end_phasic_comp_3` duration
    while (time_cyc_phasic_comp_3 > t_end_phasic_comp_3) 
        time_cyc_phasic_comp_3 = time_cyc_phasic_comp_3 - t_end_phasic_comp_3;
    end

    % If the cycle time is within the opening duration, ICC activity is present
    if (time_cyc_phasic_comp_3 <= t_open) 
        I_ICC_phasic_comp_3 = amplitude_ICC_PS;  
    else  % Otherwise, ICC activity is zero
        I_ICC_phasic_comp_3 = 0;
    end
end

%---------------------------------------------------------


%% Fundus - Stomach

% -------------------- Calcium Dynamics --------------------

% Initial intracellular calcium concentration in smooth muscle cells (SMCs) of the Fundus
Ca_i_SM_comp_1 = 0.000342013;  % (unit: mM)

% Total calcium in SMCs, incorporating tonic calcium contribution
Ca_SMC_comp_1 = Ca_i_SM_comp_1 * 1000 + Ca_fundus_tonic;  % (unit: µM)

% -------------------- Ordinary Differential Equations (ODEs) for Mechanical Activity --------------------

% These equations describe the interactions between calcium, calmodulin (CaM), myosin light chain kinase (MLCK), 
% and other biochemical components involved in smooth muscle contraction.

dy(1, :) =  Kf_Ca2CaM * Ca_SMC_comp_1^2 * CaM_comp_1 - Kr_Ca2CaM * Ca2_CaM_comp_1 ...
          - Kf_Ca2Ca4 * Ca_SMC_comp_1^2 * Ca2_CaM_comp_1 + Kr_Ca2Ca4 * Ca4_CaM_comp_1 ...
          - Kf_Ca2CaMtoCa2_CaM_MLCK * Ca2_CaM_comp_1 * MLCK_comp_1 + Kr_Ca2CaMtoCa2_CaM_MLCK * Ca2_CaM_MLCK_comp_1 ...
          + Kf_CaMBufdis * CaM_Buf_comp_1 * Ca_SMC_comp_1^2 - Kr_CaMBufdis * Ca2_CaM_comp_1 * Buf_comp_1;

dy(2, :) = -Kf_Ca2CaM * Ca_SMC_comp_1^2 * CaM_comp_1 + Kr_Ca2CaM * Ca2_CaM_comp_1 ...
          - Kf_CaMMLCK * CaM_comp_1 * MLCK_comp_1 + Kr_CaMMLCK * CaM_MLCK_comp_1 ...
          - Kf_CaMBuf * CaM_comp_1 * Buf_comp_1 + Kr_CaMBuf * CaM_Buf_comp_1;

dy(3, :) =  Kf_Ca2Ca4 * Ca_SMC_comp_1^2 * Ca2_CaM_comp_1 - Kr_Ca2Ca4 * Ca4_CaM_comp_1 ...
          - Kf_MLCKact * Ca4_CaM_comp_1 * MLCK_comp_1 + Kr_MLCKact * Ca4_CaM_MLCK_comp_1;

dy(4, :) =  Kf_CaMMLCK * CaM_comp_1 * MLCK_comp_1 - Kr_CaMMLCK * CaM_MLCK_comp_1 ...
          - Kf_Ca2CaMMLCK * Ca_SMC_comp_1^2 * CaM_MLCK_comp_1 + Kr_Ca2CaMMLCK * Ca2_CaM_MLCK_comp_1;

dy(5, :) = -Kf_CaMMLCK * CaM_comp_1 * MLCK_comp_1 + Kr_CaMMLCK * CaM_MLCK_comp_1 ...
          - Kf_Ca2CaMtoCa2_CaM_MLCK * Ca2_CaM_comp_1 * MLCK_comp_1 + Kr_Ca2CaMtoCa2_CaM_MLCK * Ca2_CaM_MLCK_comp_1 ...
          - Kf_MLCKact * Ca4_CaM_comp_1 * MLCK_comp_1 + Kr_MLCKact * Ca4_CaM_MLCK_comp_1;

dy(6, :) =  Kf_Ca2CaMtoCa2_CaM_MLCK * Ca2_CaM_comp_1 * MLCK_comp_1 - Kr_Ca2CaMtoCa2_CaM_MLCK * Ca2_CaM_MLCK_comp_1 ...
          + Kf_Ca2CaMMLCK * Ca_SMC_comp_1^2 * CaM_MLCK_comp_1 - Kr_Ca2CaMMLCK * Ca2_CaM_MLCK_comp_1 ...
          - Kf_CaCaMMLCK * Ca_SMC_comp_1^2 * Ca2_CaM_MLCK_comp_1 + Kr_CaCaMMLCK * Ca4_CaM_MLCK_comp_1;

dy(7, :) =  Kf_MLCKact * Ca4_CaM_comp_1 * MLCK_comp_1 - Kr_MLCKact * Ca4_CaM_MLCK_comp_1 ...
          + Kf_CaCaMMLCK * Ca_SMC_comp_1^2 * Ca2_CaM_MLCK_comp_1 - Kr_CaCaMMLCK * Ca4_CaM_MLCK_comp_1;

dy(8, :) =  Kf_CaMBuf * CaM_comp_1 * Buf_comp_1 - Kr_CaMBuf * CaM_Buf_comp_1 ...
          - Kf_CaMBufdis * CaM_Buf_comp_1 * Ca_SMC_comp_1^2 + Kr_CaMBufdis * Ca2_CaM_comp_1 * Buf_comp_1;

dy(9, :) = -Kf_CaMBuf * CaM_comp_1 * Buf_comp_1 + Kr_CaMBuf * CaM_Buf_comp_1 ...
          + Kf_CaMBufdis * CaM_Buf_comp_1 * Ca_SMC_comp_1^2 - Kr_CaMBufdis * Ca2_CaM_comp_1 * Buf_comp_1;

% -------------------- MLCK Activation with VIP Inhibition --------------------

% VIP reduces the amount of active Ca4_CaM_MLCK complex in the Fundus
Ca4_CaM_MLCK_comp_1 = factor_inhibitory * Ca4_CaM_MLCK_comp_1;

% -------------------- Cross-Bridge Cycling Dynamics (Hai Murphy Model)--------------------

dy(10, :) =  Vmax_MLCKact * M_comp_1 * Ca4_CaM_MLCK_comp_1 / (Km_MLCKact + M_comp_1) ...
           - HM_k3 * Mp_comp_1 + HM_k4 * AMp_comp_1 ...
           - Vmax_MLCP * Mp_comp_1 * MLCP_comp_1 / (Km_MLCP + Mp_comp_1);

dy(11, :) = -Vmax_MLCKact * M_comp_1 * Ca4_CaM_MLCK_comp_1 / (Km_MLCKact + M_comp_1) ...
           + Vmax_MLCP * Mp_comp_1 * MLCP_comp_1 / (Km_MLCP + Mp_comp_1) ...
           + HM_k7 * AM_comp_1 - HM_k8 * M_comp_1;

dy(12, :) =  Vmax_MLCKact * AM_comp_1 * Ca4_CaM_MLCK_comp_1 / (Km_MLCKact + AM_comp_1) ...
           + HM_k3 * Mp_comp_1 - HM_k4 * AMp_comp_1 ...
           - Vmax_MLCP * AMp_comp_1 * MLCP_comp_1 / (Km_MLCP + AMp_comp_1);

dy(13, :) = -Vmax_MLCKact * AM_comp_1 * Ca4_CaM_MLCK_comp_1 / (Km_MLCKact + AM_comp_1) ...
           + Vmax_MLCP * AMp_comp_1 * MLCP_comp_1 / (Km_MLCP + AMp_comp_1) ...
           - HM_k7 * AM_comp_1 + HM_k8 * M_comp_1;

% -------------------- Muscle Latch Bridge Ratio --------------------

% Compute contraction ratio based on phosphorylated and non-phosphorylated myosin
n_c_comp_1 = (AMp_comp_1 + AM_comp_1) ./ (Mp_comp_1 + M_comp_1 + AMp_comp_1 + AM_comp_1);

% Ensure numerical stability by removing complex values
n_c_comp_1 = real(n_c_comp_1);

% -------------------- Volume and Radius Adjustments --------------------

beta_vol_comp_1 = 20.689;  % Scaling factor
n_c_comp_1_max = 0.024023565;  % Maximum number of latch bridges ratio

% Compute radius adjustment factor
R_A_comp_1 = (1 - beta_vol_comp_1 * n_c_comp_1).^2 ./ (1 - beta_vol_comp_1 * n_c_comp_1_max).^2;

% Minimum volume and length ofstomach for passive stress computation due to fundus distension
volume_ini_total = 80;  
h_ini_total = 13.5;  

% Compute initial radius
radius_ini_comp_1 = sqrt(volume_ini_total / (pi * h_ini_total));

% Compute new radius based on fundus distension
R_new_comp_1 = R_A_comp_1 * radius_ini_comp_1;

% Compute new total volume based on fundus distension
V_total = pi * R_new_comp_1.^2 * h_ini_total;

% Compute passive stretch factor
const_comp_1 = 0.7 / 3.3301;
lambda_f_total = const_comp_1 * ((R_new_comp_1 / radius_ini_comp_1) - 1) + 1;

% -------------------- Compartment 1 and 2 Bridge --------------------

% Fraction of the total volume assigned to the terminal compartment
V_terminal_percentage_frac = 0.0039;

% Heights of the terminal and sphincter regions (in cm)
H_terminal = 13;   
H_sph = 13.5;      

% Compute the length of the terminal antrum compartment
L_terminal = H_sph - H_terminal;  

% Compute the radius of the terminal antrum compartment based on its toal volume of stomach by fundus distension
r_ini_comp_2 = sqrt((V_terminal_percentage_frac * V_total) / (pi * L_terminal));


%% Efferent-Afferent Mechanical and Antrum Neurosignal Processing

% -------------------- Mechanical Afferent Signaling --------------------

% Compute afferent mechanical signal based on fundus stretch factor (lambda_f_total)
f_e_afferent_mech = 42.85714286 * lambda_f_total - 42.85714286;

% Parameters for brain response to mechanical afferent signals
k_a_brain_mech = 22.77;  % Offset parameter for the brain's response function
k_b_brain_mech = 1.903;  % Scaling factor determining the steepness of the response

% -------------------- Efferent and Inhibitory Neural Signaling in the Antrum --------------------

if sympathetic == 0  
    % Compute efferent signaling to the antrum using a sigmoid function
    f_e_antrum = (0 + 0.7 * exp((f_e_afferent_mech - k_a_brain_mech) / k_b_brain_mech)) ./ ...
                 (1 + exp((f_e_afferent_mech - k_a_brain_mech) / k_b_brain_mech));
    
    % No inhibitory afferent signal in normal conditions
    f_i_antrum = 0;  
else  
    % If the sympathetic nervous system is activated, efferent signaling is suppressed
    f_e_antrum = 0;
    
    % Inhibitory input to the antrum is set to a fixed value (e.g., due to stress or sympathetic response)
    f_i_antrum = 15;  
end

%% Efferent Neuron Signaling - Antrum

% -------------------- Nitric Oxide (NO) Production --------------------

% Parameters for NO production in the antrum
a_NO_antrum = 0.667664813;
b_NO_antrum = 43.82034347;
n_NO_antrum = 1.501113355;

% Compute NO release in the antrum based on inhibitory input
NO_antrum_phasic_uM = (a_NO_antrum .* power(f_i_antrum, n_NO_antrum)) ./ ...
                       (power(b_NO_antrum, n_NO_antrum) + power(f_i_antrum, n_NO_antrum));

% Convert NO from uM to nM
NO_antrum_phasic = NO_antrum_phasic_uM * 1000;

% -------------------- cGMP Production from NO --------------------

% Parameters for cGMP production in response to NO
a_cGMP_antrum = 26.5732;
b_cGMP_antrum = 21.6111;
n_cGMP_antrum = 1.0022;

% Compute cGMP levels based on NO concentration
cGMP_antrum = (a_cGMP_antrum * power(NO_antrum_phasic, n_cGMP_antrum)) ./ ...
               (power(NO_antrum_phasic, n_cGMP_antrum) + power(b_cGMP_antrum, n_cGMP_antrum));

% Compute MLCP activation in the antrum based on cGMP levels
R_mlcp_antrum_phasic = (power(cGMP_antrum, n_mlcp)) ./ ...
                        (power(cGMP_antrum, n_mlcp) + power(K_m_mlcp, n_mlcp));

% Compute MLCP concentration in the antrum
MLCP_comp_2 = 7.5 + 2.5 * R_mlcp_antrum_phasic;

% -------------------- ICC Activity Modulation by NO --------------------

% Parameters for inhibitory modulation of ICC amplitude and frequency by NO
a_i_amplitude_ICC_antrum = 0.385470815230167;
b_i_amplitude_ICC_antrum = 0.037915896;
n_i_amplitude_ICC_antrum = 0.74976983;

a_i_cpm_antrum = 0.1228;
b_i_cpm_antrum = 3.5258;
n_i_cpm_antrum = 1.5264;

% Compute NO-mediated inhibition of ICC amplitude
phi_i_amplitude_ICC_antrum = (a_i_amplitude_ICC_antrum .* power(NO_antrum_phasic_uM, n_i_amplitude_ICC_antrum)) ./ ...
                              (power(b_i_amplitude_ICC_antrum, n_i_amplitude_ICC_antrum) + power(NO_antrum_phasic_uM, n_i_amplitude_ICC_antrum));

% Compute NO-mediated inhibition of ICC frequency
phi_i_cpm_antrum = (a_i_cpm_antrum .* power(NO_antrum_phasic_uM, n_i_cpm_antrum)) ./ ...
                    (power(b_i_cpm_antrum, n_i_cpm_antrum) + power(NO_antrum_phasic_uM, n_i_cpm_antrum));

% -------------------- Purinergic Signaling and Gap Junction Modulation --------------------

% Parameters for purinergic signaling in the antrum
a_Pur_antrum = 0.525315108793863;
b_Pur_antrum = 8.76385759972606;
n_Pur_antrum = 1.74061375397748;

% Parameters for inhibitory gap junction modulation in smooth muscle cells (SMCs)
a_i_gap_SMC_antrum = 167.5930670258801;
b_i_gap_SMC_antrum = 907.9276538870989;
n_i_gap_SMC_antrum = 0.8399709278706;

% Compute purinergic influence on the antrum based on inhibitory input
Pur_antrum = (a_Pur_antrum .* power(f_i_antrum, n_Pur_antrum)) ./ ...
             (power(b_Pur_antrum, n_Pur_antrum) + power(f_i_antrum, n_Pur_antrum));

% Compute inhibitory gap junction effect in SMCs based on purinergic signaling
phi_i_gap_SMC_antrum = (a_i_gap_SMC_antrum .* power(Pur_antrum, n_i_gap_SMC_antrum)) ./ ...
                        (power(b_i_gap_SMC_antrum, n_i_gap_SMC_antrum) + power(Pur_antrum, n_i_gap_SMC_antrum));

% -------------------- Acetylcholine (Ach) Modulation --------------------

% Parameters for Ach release in the antrum
a_ach_antrum = 156.4317944;
b_ach_antrum = 3.538564077;
n_ach_antrum = 3.993133727;

% Parameters for excitatory modulation of ICC amplitude and frequency by Ach
a_e_amplitude_ICC_antrum = 0.15;
b_e_amplitude_ICC_antrum = 0.574457529640854;
n_e_amplitude_ICC_antrum = 0.33;

a_e_cpm_antrum = 1.3;
b_e_cpm_antrum = 96;
n_e_cpm_antrum = 1;

% Compute Ach release in the antrum based on efferent input
ach_antrum = (a_ach_antrum .* power(f_e_antrum, n_ach_antrum)) ./ ...
              (power(b_ach_antrum, n_ach_antrum) + power(f_e_antrum, n_ach_antrum));

% Compute excitatory modulation of ICC amplitude by Ach
phi_e_amplitude_ICC_antrum = (a_e_amplitude_ICC_antrum .* power(ach_antrum, n_e_amplitude_ICC_antrum)) ./ ...
                              (power(b_e_amplitude_ICC_antrum, n_e_amplitude_ICC_antrum) + power(ach_antrum, n_e_amplitude_ICC_antrum));

% Compute excitatory modulation of ICC frequency by Ach
phi_e_cpm_antrum = (a_e_cpm_antrum .* power(ach_antrum, n_e_cpm_antrum)) ./ ...
                    (power(b_e_cpm_antrum, n_e_cpm_antrum) + power(ach_antrum, n_e_cpm_antrum));

% -------------------- Overall ICC Frequency and Amplitude Computation --------------------

% Compute ICC frequency (cycles per minute) incorporating excitatory and inhibitory effects
cpm_antrum = base_cpm .* (phi_e_cpm_antrum - phi_i_cpm_antrum + 1);

% Compute ICC amplitude incorporating excitatory and inhibitory effects
amplitude_ICC_antrum = base_amplitude_ICC .* (phi_e_amplitude_ICC_antrum - phi_i_amplitude_ICC_antrum + 1);

% Compute inhibitory factor from gap junction modulation in SMCs
factor_i_gap_SMC_antrum = (1 - phi_i_gap_SMC_antrum);

% -------------------- Antrum Contraction Cycle Duration --------------------

% Compute the duration of one contraction cycle in the antrum
t_end_comp_2 = 60 ./ cpm_antrum;


%% Antrum - Stomach

% Time parameters

% Start time for antrum contraction cycle (in seconds)
t_start_comp_2 = 5;

% If the cycle time is less than the start time, inhibit ICC activity
if (time_cyc_comp_2 < t_start_comp_2) 
    I_ICC_comp_2 = 0;
else 
    % Adjust the cycle time by subtracting the start time
    time_cyc_comp_2 = time_cyc_comp_2 - t_start_comp_2;

    % Restart the wave cycle after `t_end_comp_2` duration
    while (time_cyc_comp_2 > t_end_comp_2) 
        time_cyc_comp_2 = time_cyc_comp_2 - t_end_comp_2;
    end

    % If the cycle time is within the opening duration, ICC activity is present
    if (time_cyc_comp_2 <= t_open) 
        I_ICC_comp_2 = amplitude_ICC_antrum;  % ICC generates contraction signal
    else  
        I_ICC_comp_2 = 0;  % No ICC activity outside the opening duration
    end
    
end



% -------------------- ODEs for Electrical Activity --------------------

% ICC membrane potential equation based on cyclic activation
dy(14,:) = (I_ICC_comp_2 .* exp(-s_ICC.^2 .* (time_cyc_comp_2).^2) - V_m_ICC_comp_2 + V_m_rest_ICC) ./ ...
            (exp(-s_ICC.^2 .* (time_cyc_comp_2).^2));

% Compute gap junction current between ICC and SMCs
I_coup_comp_2 = G_coup .* (V_m_ICC_comp_2 - V_m_SMC_comp_2);

% SMC membrane potential equation incorporating gap junction coupling
dy(15,:) = ((factor_i_gap_SMC_antrum .* I_coup_comp_2) .* exp(-s_SMC.^2 .* (time_cyc_comp_2).^2) - V_m_SMC_comp_2 + V_m_rest_SMC) ./ ...
            (exp(-s_SMC.^2 .* (time_cyc_comp_2).^2));

% -------------------- Ionic Currents in SMCs --------------------

% Compute equilibrium potential for calcium in SMCs
E_Ca_SM_2_comp_2 = 0.5 * RToF * log(Ca_o_SM / Ca_i_SM_comp_2);

% L-type calcium channel current in SMCs
I_Ltype_SM_comp_2 = G_max_Ltype_2 * f_Ltype_SM_comp_2 * d_Ltype_SM_comp_2 * f_ca_Ltype_SM_comp_2 * (V_m_SMC_comp_2 - E_Ca_SM_2_comp_2);

% Low-voltage-activated (LVA) calcium current
E_Ca_SM_1_comp_2 = 0.5 * RToF * log(Ca_o_SM / Ca_i_SM_comp_2);
I_LVA_SM_comp_2 = G_max_LVA * f_LVA_SM_comp_2 * d_LVA_SM_comp_2 * (V_m_SMC_comp_2 - E_Ca_SM_1_comp_2);

% SR calcium release flux
J_CaSR_SM_comp_2 = J_max_CaSR * (Ca_i_SM_comp_2 * 1.0)^1.34;

% ODE for intracellular calcium concentration in SMCs
dy(16,:) = 1 * ((-I_Ltype_SM_comp_2 - I_LVA_SM_comp_2) / (2.0 * 1.0 * 1.0 * F * Vol_SM) - J_CaSR_SM_comp_2);

% -------------------- Gating Variables for Calcium Channels --------------------

% L-type calcium channel inactivation
f_inf_Ltype_SM_comp_2 = 1.0 / (1.0 + exp((V_m_SMC_comp_2 + 43.0) / 8.9));
dy(17,:) = (f_inf_Ltype_SM_comp_2 - f_Ltype_SM_comp_2) / tau_f_Ltype_SM;

% L-type calcium channel activation
d_inf_Ltype_SM_comp_2 = 1.0 / (1.0 + exp((V_m_SMC_comp_2 + 17.0) / -4.3));
dy(18,:) = (d_inf_Ltype_SM_comp_2 - d_Ltype_SM_comp_2) / tau_d_Ltype_SM;

% Calcium-dependent inactivation of L-type calcium channel
f_ca_inf_Ltype_SM_comp_2 = 1.0 - 1.0 / (1.0 + exp((Ca_i_SM_comp_2 - 0.00012401) / -0.0000131));
dy(19,:) = (f_ca_inf_Ltype_SM_comp_2 - f_ca_Ltype_SM_comp_2) / tau_f_ca_Ltype_SM;

% Low-voltage-activated (LVA) calcium channel inactivation
tau_f_LVA_SM_comp_2 = T_correction_Ca * 7.58e-3 * exp(V_m_SMC_comp_2 * 0.00817);
f_inf_LVA_SM_comp_2 = 1.0 / (1.0 + exp((V_m_SMC_comp_2 + 15.8) / 7.0));
dy(20,:) = (f_inf_LVA_SM_comp_2 - f_LVA_SM_comp_2) / tau_f_LVA_SM_comp_2;

% LVA calcium channel activation
d_inf_LVA_SM_comp_2 = 1.0 / (1.0 + exp((V_m_SMC_comp_2 + 27.5) / -10.9));
dy(21,:) = (d_inf_LVA_SM_comp_2 - d_LVA_SM_comp_2) / tau_d_LVA_SM;

% Convert intracellular calcium concentration to µM for mechanical computations
Ca_SMC_comp_2 = Ca_i_SM_comp_2 * 1000;

% -------------------- ODEs for Mechanical Activity --------------------

% Calcium-Calmodulin (CaM) and MLCK interactions regulating contraction
dy(22, :)=  Kf_Ca2CaM*Ca_SMC_comp_2*Ca_SMC_comp_2*CaM_comp_2 - Kr_Ca2CaM*Ca2_CaM_comp_2 -Kf_Ca2Ca4*Ca_SMC_comp_2*Ca_SMC_comp_2*Ca2_CaM_comp_2 + Kr_Ca2Ca4*Ca4_CaM_comp_2 -Kf_Ca2CaMtoCa2_CaM_MLCK*Ca2_CaM_comp_2*MLCK_comp_2 + Kr_Ca2CaMtoCa2_CaM_MLCK*Ca2_CaM_MLCK_comp_2 + Kf_CaMBufdis*CaM_Buf_comp_2*Ca_SMC_comp_2*Ca_SMC_comp_2 - Kr_CaMBufdis*Ca2_CaM_comp_2*Buf_comp_2;
dy(23, :)= -Kf_Ca2CaM*Ca_SMC_comp_2*Ca_SMC_comp_2*CaM_comp_2 + Kr_Ca2CaM*Ca2_CaM_comp_2 -Kf_CaMMLCK*CaM_comp_2*MLCK_comp_2 + Kr_CaMMLCK*CaM_MLCK_comp_2 -Kf_CaMBuf*CaM_comp_2*Buf_comp_2 + Kr_CaMBuf*CaM_Buf_comp_2;
dy(24, :)=  Kf_Ca2Ca4*Ca_SMC_comp_2*Ca_SMC_comp_2*Ca2_CaM_comp_2 - Kr_Ca2Ca4*Ca4_CaM_comp_2 -Kf_MLCKact*Ca4_CaM_comp_2*MLCK_comp_2 + Kr_MLCKact*Ca4_CaM_MLCK_comp_2;
dy(25, :)=  Kf_CaMMLCK*CaM_comp_2*MLCK_comp_2 - Kr_CaMMLCK*CaM_MLCK_comp_2 -Kf_Ca2CaMMLCK*Ca_SMC_comp_2*Ca_SMC_comp_2*CaM_MLCK_comp_2 + Kr_Ca2CaMMLCK*Ca2_CaM_MLCK_comp_2;
dy(26, :)= -Kf_CaMMLCK*CaM_comp_2*MLCK_comp_2 + Kr_CaMMLCK*CaM_MLCK_comp_2 -Kf_Ca2CaMtoCa2_CaM_MLCK*Ca2_CaM_comp_2*MLCK_comp_2 + Kr_Ca2CaMtoCa2_CaM_MLCK*Ca2_CaM_MLCK_comp_2 -Kf_MLCKact*Ca4_CaM_comp_2*MLCK_comp_2 + Kr_MLCKact*Ca4_CaM_MLCK_comp_2;
dy(27, :)=  Kf_Ca2CaMtoCa2_CaM_MLCK*Ca2_CaM_comp_2*MLCK_comp_2 - Kr_Ca2CaMtoCa2_CaM_MLCK*Ca2_CaM_MLCK_comp_2 + Kf_Ca2CaMMLCK*Ca_SMC_comp_2*Ca_SMC_comp_2*CaM_MLCK_comp_2 - Kr_Ca2CaMMLCK*Ca2_CaM_MLCK_comp_2 -Kf_CaCaMMLCK*Ca_SMC_comp_2*Ca_SMC_comp_2*Ca2_CaM_MLCK_comp_2 + Kr_CaCaMMLCK*Ca4_CaM_MLCK_comp_2;
dy(28, :)=  Kf_MLCKact*Ca4_CaM_comp_2*MLCK_comp_2 - Kr_MLCKact*Ca4_CaM_MLCK_comp_2 + Kf_CaCaMMLCK*Ca_SMC_comp_2*Ca_SMC_comp_2*Ca2_CaM_MLCK_comp_2 - Kr_CaCaMMLCK*Ca4_CaM_MLCK_comp_2;
dy(29, :)=  Kf_CaMBuf*CaM_comp_2*Buf_comp_2 - Kr_CaMBuf*CaM_Buf_comp_2 -Kf_CaMBufdis*CaM_Buf_comp_2*Ca_SMC_comp_2*Ca_SMC_comp_2 + Kr_CaMBufdis*Ca2_CaM_comp_2*Buf_comp_2;
dy(30, :)= -Kf_CaMBuf*CaM_comp_2*Buf_comp_2 + Kr_CaMBuf*CaM_Buf_comp_2 +Kf_CaMBufdis*CaM_Buf_comp_2*Ca_SMC_comp_2*Ca_SMC_comp_2 - Kr_CaMBufdis*Ca2_CaM_comp_2*Buf_comp_2;

% -------------------- Cross-Bridge Cycling Dynamics --------------------
%Hai Murphy Model

dy(31, :)=  Vmax_MLCKact*M_comp_2*Ca4_CaM_MLCK_comp_2/(Km_MLCKact+M_comp_2) -HM_k3*Mp_comp_2 + HM_k4*AMp_comp_2 -Vmax_MLCP*Mp_comp_2*MLCP_comp_2/(Km_MLCP+Mp_comp_2);
dy(32, :)= -Vmax_MLCKact*M_comp_2*Ca4_CaM_MLCK_comp_2/(Km_MLCKact+M_comp_2) +Vmax_MLCP*Mp_comp_2*MLCP_comp_2/(Km_MLCP+Mp_comp_2) + HM_k7*AM_comp_2 - HM_k8*M_comp_2;
dy(33, :)=  Vmax_MLCKact*AM_comp_2*Ca4_CaM_MLCK_comp_2/(Km_MLCKact+AM_comp_2) + HM_k3*Mp_comp_2 - HM_k4*AMp_comp_2 -Vmax_MLCP*AMp_comp_2*MLCP_comp_2/(Km_MLCP+AMp_comp_2);
dy(34, :)= -Vmax_MLCKact*AM_comp_2*Ca4_CaM_MLCK_comp_2/(Km_MLCKact+AM_comp_2) + Vmax_MLCP*AMp_comp_2*MLCP_comp_2/(Km_MLCP+AMp_comp_2) -HM_k7*AM_comp_2 + HM_k8*M_comp_2;


% -------------------- Mechanical Stress and Stretch --------------------

% Optimal fiber stretch properties
u_bar_fs_opt = 0.68;
x_bar_0 = 0.149;
mu_a = 72.98 * 1000;  % Convert MPa to kPa

kappa_mus = 2098.37;  % Muscle stiffness coefficient in kPa

% Compute steady-state chemical contraction effect
u_chem_bar_fs = -kappa_mus / mu_a;

% Compute contraction ratio based on phosphorylated myosin levels
n_c_comp_2 = (AMp_comp_2 + AM_comp_2) ./ (Mp_comp_2 + M_comp_2 + AMp_comp_2 + AM_comp_2);
n_c_comp_2 = real(n_c_comp_2);  % Ensure numerical stability

% Compute mechanical stretch component
u_mech_bar_fs_comp_2 = lambda_f_total - 1;

% Compute overall fiber stretch effect
u_bar_fs_comp_2 = u_mech_bar_fs_comp_2 + u_chem_bar_fs;

% Compute normalized muscle fiber length
L_bar_0_comp_2 = u_bar_fs_comp_2 - (power(u_bar_fs_comp_2, 2) / (2 * u_bar_fs_opt)) + x_bar_0;

% Compute passive element-based stress in tissue
P_a_comp_2 = mu_a * L_bar_0_comp_2 * n_c_comp_2 * (lambda_f_total - u_bar_fs_comp_2 - 1);

% Final computed tissue stress
stress_comp_2 = P_a_comp_2;

% Ensure numerical stability for stretch values
stretch_comp_2 = real(stretch_comp_2);

% -------------------- Stress-to-Stretch Relationship --------------------

% Polynomial coefficients for stretch-elasticity relation
c_5 = 1.0e+04 * 0.030597922606923;
c_4 = -0.232137790454300 * 1.0e+04;
c_3 = 0.706865564616639 * 1.0e+04;
c_2 = -1.077346195957639 * 1.0e+04;
c_1 = 0.823212261787897 * 1.0e+04;
c_0 = -0.251225167377082 * 1.0e+04;

% Compute tissue elasticity
Elast_comp_2 = c_5 .* stretch_comp_2.^5 + c_4 .* stretch_comp_2.^4 + c_3 .* stretch_comp_2.^3 + ...
               c_2 .* stretch_comp_2.^2 + c_1 .* stretch_comp_2.^1 + c_0;  

% Compute mechanical damping (dashpot effect)
dashpot_1_comp_2 = S_4_comp_2_2_3_4 .* exp(S_5_comp_2_2_3_4 .* (stress_comp_2 - Elast_comp_2 * S_3_comp_2_2_3_4));
dashpot_2_comp_2 = dashpot_1_comp_2 .* tanh(S_6_comp_2_2_3_4 .* power(stretch_comp_2 - 0.98, 2));

% Compute stretch velocity
dy(35,:) = (1 / dashpot_2_comp_2) * (stress_comp_2 - Elast_comp_2);

% Compute compartment radius
[r_fin_comp_2, r_d_comp_2] = compartment_radius(stretch_comp_2, L_t_comp_2, r_ini_comp_2);



%% Pyloric Sphincter - Stomach

% -------------------- Calcium and MLCP Dynamics --------------------

% Initial intracellular calcium concentration in smooth muscle cells (SMCs) of the pyloric sphincter (tonic component)
Ca_i_SM_tonic_comp_3 = 1.817624401818595e-04;  % (unit: mM)

% Compute total calcium concentration in SMCs, incorporating Ach-mediated calcium release
Ca_SMC_tonic_comp_3 = Ca_i_SM_tonic_comp_3 * 1000 + Ca_SMC_PS_ach_tonic;  % (unit: µM)

% Compute MLCP activation in the pyloric sphincter based on non-adrenergic, non-cholinergic (NANC) signaling
MLCP_PS_tonic = 7.5 + 2.5 * factor_NANC_PS_tonic;

% -------------------- ODEs for Mechanical Activity --------------------

% These equations describe calcium binding to calmodulin (CaM), myosin light chain kinase (MLCK) activation,
% and other biochemical processes regulating smooth muscle contraction in the pyloric sphincter.


dy(36, :)=  Kf_Ca2CaM*Ca_SMC_tonic_comp_3*Ca_SMC_tonic_comp_3*CaM_tonic_comp_3 - Kr_Ca2CaM*Ca2_CaM_tonic_comp_3 -Kf_Ca2Ca4*Ca_SMC_tonic_comp_3*Ca_SMC_tonic_comp_3*Ca2_CaM_tonic_comp_3 + Kr_Ca2Ca4*Ca4_CaM_tonic_comp_3 -Kf_Ca2CaMtoCa2_CaM_MLCK*Ca2_CaM_tonic_comp_3*MLCK_tonic_comp_3 + Kr_Ca2CaMtoCa2_CaM_MLCK*Ca2_CaM_MLCK_tonic_comp_3 + Kf_CaMBufdis*CaM_Buf_tonic_comp_3*Ca_SMC_tonic_comp_3*Ca_SMC_tonic_comp_3 - Kr_CaMBufdis*Ca2_CaM_tonic_comp_3*Buf_tonic_comp_3;
dy(37, :)= -Kf_Ca2CaM*Ca_SMC_tonic_comp_3*Ca_SMC_tonic_comp_3*CaM_tonic_comp_3 + Kr_Ca2CaM*Ca2_CaM_tonic_comp_3 -Kf_CaMMLCK*CaM_tonic_comp_3*MLCK_tonic_comp_3 + Kr_CaMMLCK*CaM_MLCK_tonic_comp_3 -Kf_CaMBuf*CaM_tonic_comp_3*Buf_tonic_comp_3 + Kr_CaMBuf*CaM_Buf_tonic_comp_3;
dy(38, :)=  Kf_Ca2Ca4*Ca_SMC_tonic_comp_3*Ca_SMC_tonic_comp_3*Ca2_CaM_tonic_comp_3 - Kr_Ca2Ca4*Ca4_CaM_tonic_comp_3 -Kf_MLCKact*Ca4_CaM_tonic_comp_3*MLCK_tonic_comp_3 + Kr_MLCKact*Ca4_CaM_MLCK_tonic_comp_3;
dy(39, :)=  Kf_CaMMLCK*CaM_tonic_comp_3*MLCK_tonic_comp_3 - Kr_CaMMLCK*CaM_MLCK_tonic_comp_3 -Kf_Ca2CaMMLCK*Ca_SMC_tonic_comp_3*Ca_SMC_tonic_comp_3*CaM_MLCK_tonic_comp_3 + Kr_Ca2CaMMLCK*Ca2_CaM_MLCK_tonic_comp_3;
dy(40, :)= -Kf_CaMMLCK*CaM_tonic_comp_3*MLCK_tonic_comp_3 + Kr_CaMMLCK*CaM_MLCK_tonic_comp_3 -Kf_Ca2CaMtoCa2_CaM_MLCK*Ca2_CaM_tonic_comp_3*MLCK_tonic_comp_3 + Kr_Ca2CaMtoCa2_CaM_MLCK*Ca2_CaM_MLCK_tonic_comp_3 -Kf_MLCKact*Ca4_CaM_tonic_comp_3*MLCK_tonic_comp_3 + Kr_MLCKact*Ca4_CaM_MLCK_tonic_comp_3;
dy(41, :)=  Kf_Ca2CaMtoCa2_CaM_MLCK*Ca2_CaM_tonic_comp_3*MLCK_tonic_comp_3 - Kr_Ca2CaMtoCa2_CaM_MLCK*Ca2_CaM_MLCK_tonic_comp_3 + Kf_Ca2CaMMLCK*Ca_SMC_tonic_comp_3*Ca_SMC_tonic_comp_3*CaM_MLCK_tonic_comp_3 - Kr_Ca2CaMMLCK*Ca2_CaM_MLCK_tonic_comp_3 -Kf_CaCaMMLCK*Ca_SMC_tonic_comp_3*Ca_SMC_tonic_comp_3*Ca2_CaM_MLCK_tonic_comp_3 + Kr_CaCaMMLCK*Ca4_CaM_MLCK_tonic_comp_3;
dy(42, :)=  Kf_MLCKact*Ca4_CaM_tonic_comp_3*MLCK_tonic_comp_3 - Kr_MLCKact*Ca4_CaM_MLCK_tonic_comp_3 + Kf_CaCaMMLCK*Ca_SMC_tonic_comp_3*Ca_SMC_tonic_comp_3*Ca2_CaM_MLCK_tonic_comp_3 - Kr_CaCaMMLCK*Ca4_CaM_MLCK_tonic_comp_3;
dy(43, :)=  Kf_CaMBuf*CaM_tonic_comp_3*Buf_tonic_comp_3 - Kr_CaMBuf*CaM_Buf_tonic_comp_3 -Kf_CaMBufdis*CaM_Buf_tonic_comp_3*Ca_SMC_tonic_comp_3*Ca_SMC_tonic_comp_3 + Kr_CaMBufdis*Ca2_CaM_tonic_comp_3*Buf_tonic_comp_3;
dy(44, :)= -Kf_CaMBuf*CaM_tonic_comp_3*Buf_tonic_comp_3 + Kr_CaMBuf*CaM_Buf_tonic_comp_3 +Kf_CaMBufdis*CaM_Buf_tonic_comp_3*Ca_SMC_tonic_comp_3*Ca_SMC_tonic_comp_3 - Kr_CaMBufdis*Ca2_CaM_tonic_comp_3*Buf_tonic_comp_3;


% -------------------- MLCK Activation with Inhibitory Modulation --------------------

% Compute inhibition factor based on NANC signaling
factor_inh_PS_tonic = 1 - factor_NANC_PS_tonic;

% Compute the effective Ca4-CaM-MLCK complex level, incorporating inhibitory effects
Ca4_CaM_MLCK_tonic_comp_3 = factor_inh_PS_tonic * Ca4_CaM_MLCK_tonic_comp_3;

% -------------------- Latch-Bridge Cycling Dynamics --------------------

% Myosin phosphorylation and MLCK regulation in the tonic component of the pyloric sphincter

%Hai Murphy Model
dy(45, :)=  Vmax_MLCKact*M_tonic_comp_3*Ca4_CaM_MLCK_tonic_comp_3/(Km_MLCKact+M_tonic_comp_3) -HM_k3*Mp_tonic_comp_3 + HM_k4*AMp_tonic_comp_3 -Vmax_MLCP*Mp_tonic_comp_3*MLCP_PS_tonic/(Km_MLCP+Mp_tonic_comp_3);
dy(46, :)= -Vmax_MLCKact*M_tonic_comp_3*Ca4_CaM_MLCK_tonic_comp_3/(Km_MLCKact+M_tonic_comp_3) +Vmax_MLCP*Mp_tonic_comp_3*MLCP_PS_tonic/(Km_MLCP+Mp_tonic_comp_3) + HM_k7*AM_tonic_comp_3 - HM_k8*M_tonic_comp_3;
dy(47, :)=  Vmax_MLCKact*AM_tonic_comp_3*Ca4_CaM_MLCK_tonic_comp_3/(Km_MLCKact+AM_tonic_comp_3) + HM_k3*Mp_tonic_comp_3 - HM_k4*AMp_tonic_comp_3 -Vmax_MLCP*AMp_tonic_comp_3*MLCP_PS_tonic/(Km_MLCP+AMp_tonic_comp_3);
dy(48, :)= -Vmax_MLCKact*AM_tonic_comp_3*Ca4_CaM_MLCK_tonic_comp_3/(Km_MLCKact+AM_tonic_comp_3) + Vmax_MLCP*AMp_tonic_comp_3*MLCP_PS_tonic/(Km_MLCP+AMp_tonic_comp_3) -HM_k7*AM_tonic_comp_3 + HM_k8*M_tonic_comp_3;


% -------------------- Number of latch bridges ratio --------------------

% Compute contraction ratio based on phosphorylated and non-phosphorylated myosin levels
n_c_tonic_comp_3 = (AMp_tonic_comp_3 + AM_tonic_comp_3) ./ ...
                    (Mp_tonic_comp_3 + M_tonic_comp_3 + AMp_tonic_comp_3 + AM_tonic_comp_3);

% Ensure numerical stability by removing complex values
n_c_tonic_comp_3 = real(n_c_tonic_comp_3);

% -------------------- Radius and Stretch Adjustments --------------------

% Scaling factor for contraction-volume relationship
beta_vol_tonic_PS = 236.316320000;

% Maximum contraction ratio
n_c_tonic_comp_3_max = 0.002559025607741;

% Compute radius adjustment factor
R_A_tonic_PS = (1 - beta_vol_tonic_PS * n_c_tonic_comp_3).^2 ./ ...
               (1 - beta_vol_tonic_PS * n_c_tonic_comp_3_max).^2;

% Initial radius of the pyloric sphincter in the tonic state
radius_ini_tonic_PS = 0.075;  % (unit: cm)

% Compute new radius based on PS distension
R_new_tonic_PS = R_A_tonic_PS .* radius_ini_tonic_PS;

% Assign the tonic radius to the phasic compartment
r_ini_phasic_comp_3 = R_new_tonic_PS;

% Prevent radius from exceeding physiological limits
if r_ini_phasic_comp_3 <= 0
    r_ini_phasic_comp_3 = 0;
end

if r_ini_phasic_comp_3 >= 0.47
    r_ini_phasic_comp_3 = 0.47;
end

% -------------------- Stretch Computation --------------------

% Scaling constant for stretch calculation
const_tonic_PS = 0.1296;

% Compute the passive stretch factor based on radius changes
lambda_f_PS = const_tonic_PS .* ((R_new_tonic_PS ./ radius_ini_tonic_PS) - 1) + 1;


%% Phasic Section - Pyloric Sphincter

% -------------------- ODEs for Electrical Activity --------------------

% ICC membrane potential equation based on cyclic activation in the phasic compartment
dy(49,:) = (I_ICC_phasic_comp_3 .* exp(-s_ICC.^2 .* (time_cyc_phasic_comp_3).^2) - V_m_ICC_phasic_comp_3 + V_m_rest_ICC) ./ ...
            (exp(-s_ICC.^2 .* (time_cyc_phasic_comp_3).^2));

% Compute gap junction current between ICC and SMCs
I_coup_phasic_comp_3 = G_coup .* (V_m_ICC_phasic_comp_3 - V_m_SMC_phasic_comp_3);

% SMC membrane potential equation incorporating gap junction coupling
dy(50,:) = ((factor_i_gap_SMC_PS .* I_coup_phasic_comp_3) .* exp(-s_SMC.^2 .* (time_cyc_phasic_comp_3).^2) - V_m_SMC_phasic_comp_3 + V_m_rest_SMC) ./ ...
            (exp(-s_SMC.^2 .* (time_cyc_phasic_comp_3).^2));

% -------------------- Ionic Currents in SMCs --------------------

% Compute equilibrium potential for calcium in SMCs
E_Ca_SM_2_phasic_comp_3 = 0.5 * RToF * log(Ca_o_SM / Ca_i_SM_phasic_comp_3);

% L-type calcium channel current in SMCs
I_Ltype_SM_phasic_comp_3 = G_max_Ltype_2 * f_Ltype_SM_phasic_comp_3 * d_Ltype_SM_phasic_comp_3 * f_ca_Ltype_SM_phasic_comp_3 * ...
                            (V_m_SMC_phasic_comp_3 - E_Ca_SM_2_phasic_comp_3);

% Low-voltage-activated (LVA) calcium current
E_Ca_SM_1_phasic_comp_3 = 0.5 * RToF * log(Ca_o_SM / Ca_i_SM_phasic_comp_3);
I_LVA_SM_phasic_comp_3 = G_max_LVA * f_LVA_SM_phasic_comp_3 * d_LVA_SM_phasic_comp_3 * (V_m_SMC_phasic_comp_3 - E_Ca_SM_1_phasic_comp_3);

% SR calcium release flux
J_CaSR_SM_phasic_comp_3 = J_max_CaSR * (Ca_i_SM_phasic_comp_3 * 1.0)^1.34;

% ODE for intracellular calcium concentration in SMCs
dy(51,:) = 1 * ((-I_Ltype_SM_phasic_comp_3 - I_LVA_SM_phasic_comp_3) / (2.0 * 1.0 * 1.0 * F * Vol_SM) - J_CaSR_SM_phasic_comp_3);

% -------------------- Gating Variables for Calcium Channels --------------------

% L-type calcium channel inactivation
f_inf_Ltype_SM_phasic_comp_3 = 1.0 / (1.0 + exp((V_m_SMC_phasic_comp_3 + 43.0) / 8.9));
dy(52,:) = (f_inf_Ltype_SM_phasic_comp_3 - f_Ltype_SM_phasic_comp_3) / tau_f_Ltype_SM;

% L-type calcium channel activation
d_inf_Ltype_SM_phasic_comp_3 = 1.0 / (1.0 + exp((V_m_SMC_phasic_comp_3 + 17.0) / -4.3));
dy(53,:) = (d_inf_Ltype_SM_phasic_comp_3 - d_Ltype_SM_phasic_comp_3) / tau_d_Ltype_SM;

% Calcium-dependent inactivation of L-type calcium channel
f_ca_inf_Ltype_SM_phasic_comp_3 = 1.0 - 1.0 / (1.0 + exp((Ca_i_SM_phasic_comp_3 - 0.00012401) / -0.0000131));
dy(54,:) = (f_ca_inf_Ltype_SM_phasic_comp_3 - f_ca_Ltype_SM_phasic_comp_3) / tau_f_ca_Ltype_SM;

% Low-voltage-activated (LVA) calcium channel inactivation
tau_f_LVA_SM_phasic_comp_3 = T_correction_Ca * 7.58e-3 * exp(V_m_SMC_phasic_comp_3 * 0.00817);
f_inf_LVA_SM_phasic_comp_3 = 1.0 / (1.0 + exp((V_m_SMC_phasic_comp_3 + 15.8) / 7.0));
dy(55,:) = (f_inf_LVA_SM_phasic_comp_3 - f_LVA_SM_phasic_comp_3) / tau_f_LVA_SM_phasic_comp_3;

% LVA calcium channel activation
d_inf_LVA_SM_phasic_comp_3 = 1.0 / (1.0 + exp((V_m_SMC_phasic_comp_3 + 27.5) / -10.9));
dy(56,:) = (d_inf_LVA_SM_phasic_comp_3 - d_LVA_SM_phasic_comp_3) / tau_d_LVA_SM;

% Convert intracellular calcium concentration to µM for mechanical computations
Ca_SMC_phasic_comp_3 = Ca_i_SM_phasic_comp_3 * 1000;

% -------------------- ODEs for Mechanical Activity Ca - > MLCK active calculation --------------------

dy(57, :)=  Kf_Ca2CaM*Ca_SMC_phasic_comp_3*Ca_SMC_phasic_comp_3*CaM_phasic_comp_3 - Kr_Ca2CaM*Ca2_CaM_phasic_comp_3 -Kf_Ca2Ca4*Ca_SMC_phasic_comp_3*Ca_SMC_phasic_comp_3*Ca2_CaM_phasic_comp_3 + Kr_Ca2Ca4*Ca4_CaM_phasic_comp_3 -Kf_Ca2CaMtoCa2_CaM_MLCK*Ca2_CaM_phasic_comp_3*MLCK_phasic_comp_3 + Kr_Ca2CaMtoCa2_CaM_MLCK*Ca2_CaM_MLCK_phasic_comp_3 + Kf_CaMBufdis*CaM_Buf_phasic_comp_3*Ca_SMC_phasic_comp_3*Ca_SMC_phasic_comp_3 - Kr_CaMBufdis*Ca2_CaM_phasic_comp_3*Buf_phasic_comp_3;
dy(58, :)= -Kf_Ca2CaM*Ca_SMC_phasic_comp_3*Ca_SMC_phasic_comp_3*CaM_phasic_comp_3 + Kr_Ca2CaM*Ca2_CaM_phasic_comp_3 -Kf_CaMMLCK*CaM_phasic_comp_3*MLCK_phasic_comp_3 + Kr_CaMMLCK*CaM_MLCK_phasic_comp_3 -Kf_CaMBuf*CaM_phasic_comp_3*Buf_phasic_comp_3 + Kr_CaMBuf*CaM_Buf_phasic_comp_3;
dy(59, :)=  Kf_Ca2Ca4*Ca_SMC_phasic_comp_3*Ca_SMC_phasic_comp_3*Ca2_CaM_phasic_comp_3 - Kr_Ca2Ca4*Ca4_CaM_phasic_comp_3 -Kf_MLCKact*Ca4_CaM_phasic_comp_3*MLCK_phasic_comp_3 + Kr_MLCKact*Ca4_CaM_MLCK_phasic_comp_3;
dy(60, :)=  Kf_CaMMLCK*CaM_phasic_comp_3*MLCK_phasic_comp_3 - Kr_CaMMLCK*CaM_MLCK_phasic_comp_3 -Kf_Ca2CaMMLCK*Ca_SMC_phasic_comp_3*Ca_SMC_phasic_comp_3*CaM_MLCK_phasic_comp_3 + Kr_Ca2CaMMLCK*Ca2_CaM_MLCK_phasic_comp_3;
dy(61, :)= -Kf_CaMMLCK*CaM_phasic_comp_3*MLCK_phasic_comp_3 + Kr_CaMMLCK*CaM_MLCK_phasic_comp_3 -Kf_Ca2CaMtoCa2_CaM_MLCK*Ca2_CaM_phasic_comp_3*MLCK_phasic_comp_3 + Kr_Ca2CaMtoCa2_CaM_MLCK*Ca2_CaM_MLCK_phasic_comp_3 -Kf_MLCKact*Ca4_CaM_phasic_comp_3*MLCK_phasic_comp_3 + Kr_MLCKact*Ca4_CaM_MLCK_phasic_comp_3;
dy(62, :)=  Kf_Ca2CaMtoCa2_CaM_MLCK*Ca2_CaM_phasic_comp_3*MLCK_phasic_comp_3 - Kr_Ca2CaMtoCa2_CaM_MLCK*Ca2_CaM_MLCK_phasic_comp_3 + Kf_Ca2CaMMLCK*Ca_SMC_phasic_comp_3*Ca_SMC_phasic_comp_3*CaM_MLCK_phasic_comp_3 - Kr_Ca2CaMMLCK*Ca2_CaM_MLCK_phasic_comp_3 -Kf_CaCaMMLCK*Ca_SMC_phasic_comp_3*Ca_SMC_phasic_comp_3*Ca2_CaM_MLCK_phasic_comp_3 + Kr_CaCaMMLCK*Ca4_CaM_MLCK_phasic_comp_3;
dy(63, :)=  Kf_MLCKact*Ca4_CaM_phasic_comp_3*MLCK_phasic_comp_3 - Kr_MLCKact*Ca4_CaM_MLCK_phasic_comp_3 + Kf_CaCaMMLCK*Ca_SMC_phasic_comp_3*Ca_SMC_phasic_comp_3*Ca2_CaM_MLCK_phasic_comp_3 - Kr_CaCaMMLCK*Ca4_CaM_MLCK_phasic_comp_3;
dy(64, :)=  Kf_CaMBuf*CaM_phasic_comp_3*Buf_phasic_comp_3 - Kr_CaMBuf*CaM_Buf_phasic_comp_3 -Kf_CaMBufdis*CaM_Buf_phasic_comp_3*Ca_SMC_phasic_comp_3*Ca_SMC_phasic_comp_3 + Kr_CaMBufdis*Ca2_CaM_phasic_comp_3*Buf_phasic_comp_3;
dy(65, :)= -Kf_CaMBuf*CaM_phasic_comp_3*Buf_phasic_comp_3 + Kr_CaMBuf*CaM_Buf_phasic_comp_3 +Kf_CaMBufdis*CaM_Buf_phasic_comp_3*Ca_SMC_phasic_comp_3*Ca_SMC_phasic_comp_3 - Kr_CaMBufdis*Ca2_CaM_phasic_comp_3*Buf_phasic_comp_3;


% -------------------- Cross-Bridge Cycling Dynamics --------------------

% Myosin phosphorylation and MLCK regulation in the phasic component of the pyloric sphincter

%Hai Murphy model
dy(66, :)=  Vmax_MLCKact*M_phasic_comp_3*Ca4_CaM_MLCK_phasic_comp_3/(Km_MLCKact+M_phasic_comp_3) -HM_k3*Mp_phasic_comp_3 + HM_k4*AMp_phasic_comp_3 -Vmax_MLCP*Mp_phasic_comp_3*MLCP_phasic_comp_3/(Km_MLCP+Mp_phasic_comp_3);
dy(67, :)= -Vmax_MLCKact*M_phasic_comp_3*Ca4_CaM_MLCK_phasic_comp_3/(Km_MLCKact+M_phasic_comp_3) +Vmax_MLCP*Mp_phasic_comp_3*MLCP_phasic_comp_3/(Km_MLCP+Mp_phasic_comp_3) + HM_k7*AM_phasic_comp_3 - HM_k8*M_phasic_comp_3;
dy(68, :)=  Vmax_MLCKact*AM_phasic_comp_3*Ca4_CaM_MLCK_phasic_comp_3/(Km_MLCKact+AM_phasic_comp_3) + HM_k3*Mp_phasic_comp_3 - HM_k4*AMp_phasic_comp_3 -Vmax_MLCP*AMp_phasic_comp_3*MLCP_phasic_comp_3/(Km_MLCP+AMp_phasic_comp_3);
dy(69, :)= -Vmax_MLCKact*AM_phasic_comp_3*Ca4_CaM_MLCK_phasic_comp_3/(Km_MLCKact+AM_phasic_comp_3) + Vmax_MLCP*AMp_phasic_comp_3*MLCP_phasic_comp_3/(Km_MLCP+AMp_phasic_comp_3) -HM_k7*AM_phasic_comp_3 + HM_k8*M_phasic_comp_3;

% -------------------- Muscle Fiber Stretch and Stress --------------------

% Optimal fiber stretch properties
u_bar_fs_opt = 0.68;  % Optimal fiber stretch ratio
x_bar_0 = 0.149;      % Offset parameter for stretch computation
mu_a = 72.98 * 1000;  % Convert MPa to kPa (material stiffness)

% Muscle stiffness coefficient
kappa_mus = 2098.37;  % kPa

% Compute steady-state chemical contraction effect
u_chem_bar_fs = -kappa_mus / mu_a;

% -------------------- Contraction Ratio Computation --------------------

% Compute contraction ratio based on phosphorylated myosin levels
n_c_phasic_comp_3 = (AMp_phasic_comp_3 + AM_phasic_comp_3) ./ ...
                     (Mp_phasic_comp_3 + M_phasic_comp_3 + AMp_phasic_comp_3 + AM_phasic_comp_3);

% Ensure numerical stability by removing complex values
n_c_phasic_comp_3 = real(n_c_phasic_comp_3);

% Compute mechanical stretch component
u_mech_bar_fs_phasic_comp_3 = lambda_f_PS - 1;

% Compute overall fiber stretch effect
u_bar_fs_phasic_comp_3 = u_mech_bar_fs_phasic_comp_3 + u_chem_bar_fs;

% Compute normalized muscle fiber length
L_bar_0_phasic_comp_3 = u_bar_fs_phasic_comp_3 - (power(u_bar_fs_phasic_comp_3, 2) / (2 * u_bar_fs_opt)) + x_bar_0;

% -------------------- Stress Computation --------------------

% Compute passive element-based stress in tissue
P_a_phasic_comp_3 = mu_a * L_bar_0_phasic_comp_3 * n_c_phasic_comp_3 * (lambda_f_PS - u_bar_fs_phasic_comp_3 - 1);

% Prevent negative stress values (ensuring non-compressive forces)
if P_a_phasic_comp_3 <= 0
    P_a_phasic_comp_3 = 0;
end

% Assign final computed tissue stress
stress_phasic_comp_3 = P_a_phasic_comp_3;

% Ensure numerical stability for stretch values
stretch_phasic_comp_3 = real(stretch_phasic_comp_3);

% -------------------- Stress-to-Stretch Relationship --------------------

% Compute tissue elasticity response to stretch
Elast_phasic_comp_3 = c_5 .* stretch_phasic_comp_3.^5 + c_4 .* stretch_phasic_comp_3.^4 + ...
                      c_3 .* stretch_phasic_comp_3.^3 + c_2 .* stretch_phasic_comp_3.^2 + ...
                      c_1 .* stretch_phasic_comp_3 + c_0;  

% -------------------- Dashpot Model (Viscoelastic Damping) --------------------

% Compute first damping coefficient based on tissue stress
dashpot_1_phasic_comp_3 = S_4_phasic_comp_3_2_3_4 .* exp(S_5_phasic_comp_3_2_3_4 .* ...
                                (stress_phasic_comp_3 - Elast_phasic_comp_3 * S_3_phasic_comp_3_2_3_4));

% Compute second damping coefficient based on stretch magnitude
dashpot_2_phasic_comp_3 = dashpot_1_phasic_comp_3 .* tanh(S_6_phasic_comp_3_2_3_4 .* ...
                                power(stretch_phasic_comp_3 - 0.98, 2));

% Compute stretch velocity using viscoelastic properties
dy(70, :) = (1 / dashpot_2_phasic_comp_3) * (stress_phasic_comp_3 - Elast_phasic_comp_3);

% -------------------- Compartment Radius Computation --------------------

% Compute new compartment radius based on stretch
[r_fin_phasic_comp_3, r_d_phasic_comp_3] = compartment_radius(stretch_phasic_comp_3, L_t_phasic_comp_3, r_ini_phasic_comp_3);

% Prevent non-physiological radius values
if r_fin_phasic_comp_3 <= 0.0001
    r_fin_phasic_comp_3 = 0;
end

% -------------------- Flow Rate Computation --------------------

% Maximum flow rate through the compartment (converted to per second)
Q_max = 70 / 60;  % Convert from ml/min to ml/sec

% Compute flow rate based on radius using Poiseuille's approximation (proportional to radius squared)
Q_flow = Q_max .* power(r_fin_phasic_comp_3 / 0.467644, 2);

% Set flow rate to zero for the initial phase of simulation
if time <= 5
    Q_flow = 0;
end

% -------------------- Flow Rate Change Over Time --------------------

% Compute the time derivative of fluid flow within the compartment.
% This represents the rate of volume change due to fluid movement.

dy(71, :) = -Q_flow;  % Flow rate reduction over time



end