%% Atmosphere Function
function [Z, Z_L, Z_U, T, P, rho, c, g, mu, nu, k, n, n_sum] = atmo(alt,division,units)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Program:    1976 Standard Atmosphere Calculator[0-1000 km]
%   Author:     Brent Lewis(RocketLion@gmail.com)
%               University of Colorado-Boulder
%   History:    Original-1/10/2007
%               Revision-1/12/2007-Corrected for changes in Matlab versions
%               for backward compatability-Many thanks to Rich
%               Rieber(rrieber@gmail.com)
%   Input:      alt:        Final Geometric Altitude[km]
%               division:   Reporting points for output arrays[km]
%                           (.01 km & Divisible by .01 km)
%               units:      1-[Metric]
%                           2-{English}
%   Default:    Values used if no input
%               alt:        1000 km
%               division:   1 km
%               units:      Metric
%   Output:     Each value has a specific region that it is valid in with this model
%               and is only printed out in that region
%               Z:          Total Reporting Altitudes[0<=alt<=1000 km][km]{ft}
%               Z_L:        Lower Atmosphere Reporting Altitudes[0<=alt<=86 km][km]{ft}
%               Z_U:        Upper Atmosphere Reporting Altitudes[86<=alt<=1000 km][km]{ft}
%               T:          Temperature array[0<=alt<=1000 km][K]{R}
%               P:          Pressure array[0<=alt<=1000 km][Pa]{in_Hg}
%               rho:        Density array[0<=alt<=1000 km][kg/m^3]{lb/ft^3}
%               c:          Speed of sound array[0<=alt<=86 km][m/s]{ft/s}
%               g:          Gravity array[0<=alt<=1000 km][m/s^2]{ft/s^2}
%               mu:         Dynamic Viscosity array[0<=alt<=86 km][N*s/m^2]{lb/(ft*s)}
%               nu:         Kinematic Viscosity array[0<=alt<=86 km][m^2/s]{ft^2/s}
%               k:          Coefficient of Thermal Conductivity
%                           array[0<=alt<=86 km][W/(m*K)]{BTU/(ft*s*R)}
%               n:          Number Density of individual gases
%                           (N2 O O2 Ar He H)[86km<=alt<=1000km][1/m^3]{1/ft^3}
%               n_sum:      Number Density of total gases
%                           [86km<=alt<=1000km][1/m^3]{1/ft^3}
%   Acknowledgements:       1976 U.S. Standard Atmosphere
%                           Prof. Adam Norris-Numerical Analysis Class
%                           Steven S. Pietrobon USSA1976 Program
%   Notes:                  Program uses a 5-point Simpsons Rule in 10
%                           meter increments.  Results DO vary by less 1%
%                           compared to tabulated values and is probably
%                           caused by different integration techniques
%   Examples:               atmo() will compute the full atmosphere in 1 km
%                           increments and output in Metric Units
%                           atmo(10) will compute the atmosphere between 0
%                           and 10 km in 1 km increments and output in
%                           Metric Units
%                           atmo(20,.1,2) will compute the atmosphere
%                           between 0 and 20 km in 100 m increments and
%                           output in English Units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

if nargin == 0
    alt = 1000;
    division = 1;
    units = 1;
elseif nargin == 1
    division = 1;
    units = 1;
elseif nargin == 2
    units = 1;
end

%   Error Reporting
if nargin > 3
    error('Too many inputs')
elseif mod(division,.01) ~= 0
    error('Divisions must be multiples of .01 km')
elseif units ~= 1 && units ~= 2
    error('Units Choice Invalid[1-Metric,2-English]')
elseif alt<0 || alt>1000
    error('Program only valid for 0<altitudes<1000 km')
end

%   Matrix Pre-allocation
if alt <= 86
    Z_L = (0:division:alt);
    Z_U = [];
    n = [];
else
    Z_L = (0:division:86)';
    Z_U = (86:division:alt)';
    if mod(86,division) ~= 0
        Z_L = [Z_L; 86];
    end
    if mod(alt-86,division) ~= 0
        Z_U = [Z_U; alt];
    end
end
T_L = zeros(size(Z_L));
T_M_L = T_L;
T_U = zeros(size(Z_U));

%   Conversion Factor Used in 80<alt<86 km
Z_M = 80:.5:86;
M_M_0 = [1 .999996 .999989 .999971 .999941 .999909 ...
    .999870 .999829 .999786 .999741 .999694 .999641 .999579];

%   Constants
M_0 = 28.9644;
M_i = [28.0134; 15.9994; 31.9988; 39.948; 4.0026; 1.00797];
beta = 1.458e-6;
gamma = 1.4;
g_0 = 9.80665;
R = 8.31432e3;
r_E = 6.356766e3;
S = 110.4;
N_A = 6.022169e26;

%   Temperature
for i = 1 : length(Z_L)
    T_L(i,1) = atmo_temp(Z_L(i));
    T_M_L(i,1) = T_L(i,1);
    if Z_L(i) > 80 && Z_L(i) < 86
        T_L(i,1) = T_L(i)*interp1(Z_M,M_M_0,Z_L(i));
    end
end
for i = 1 : length(Z_U)
    T_U(i,1) = atmo_temp(Z_U(i));
end

%   Number Density
if alt > 86
    n = atmo_compo(alt,division);
    n_sum = sum(n,2);
else
    n = [];
    n_sum = [];
end

%   Pressure
P_L = atmo_p(Z_L);
P_U = atmo_p(Z_U,T_U,n_sum);

%   Density
rho_L = M_0*P_L./(R*T_M_L);
if ~isempty(P_U)
    rho_U = n*M_i/N_A;
else
    rho_U = [];
end

%   Speed of Sound
c = sqrt(gamma*R*T_M_L/M_0);
%   Dynamic Viscosity
mu = beta*T_L.^1.5./(T_L+S);
%   Kinematic Viscosity
nu = mu./rho_L;
%   Thermal Conductivity Coefficient
k = 2.64638e-3*T_L.^1.5./(T_L+245*10.^(-12./T_L));

%   Combine Models
T = [T_L(1:end-1*double(~isempty(T_U)));T_U];
P = [P_L(1:end-1*double(~isempty(T_U)));P_U];
rho = [rho_L(1:end-1*double(~isempty(T_U)));rho_U];
Z = [Z_L(1:end-1*double(~isempty(T_U)));Z_U];

%   Gravity
g = g_0*(r_E./(r_E+Z)).^2;

if units == 2
    unit_c = [3.048e-1 3.048e-1 3.048e-1 5/9 0.0001450377 1.6018463e1...
        3.048e-1 3.048e-1 1.488163944 9.290304e-2 6.226477504e-3...
        3.531466672e2 3.531466672e2];
    Z = Z/unit_c(1);
    Z_L = Z_L/unit_c(2);
    Z_U = Z_U/unit_c(3);
    T = T/unit_c(4);
    P = P/unit_c(5);
    rho = rho/unit_c(6);
    c = c/unit_c(7); 
    g = g/unit_c(8);
    mu = mu/unit_c(9);
    nu = nu/unit_c(10); 
    k = n/unit_c(11);
    n_sum = n_sum/unit_c(12);
end
end

%% Composition Function
function n_i_array = atmo_compo(alt,division)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Program:    High Altitude Atmospheric Composition Calculation
%   Author:     Brent Lewis(RocketLion@gmail.com)
%               University of Colorado-Boulder
%   History:    Original-1/10/2007
%               Revision-1/12/2007-Corrected for changes in Matlab versions
%               for backward compatability-Many thanks to Rich
%               Rieber(rrieber@gmail.com)
%   Input:      alt:        Geometric Altitude of desired altitude[scalar][km] 
%               division:   Desired output altitudes
%   Output:     n_i_array:  Array of compositions of [N2 O O2 Ar He H] at
%                           desired reporting altitudes using equations
%                           from 1976 Standard Atmosphere
%   Note:       Only Valid Between 86 km and 1000 km
%               Division must be a multiple of 10 m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z_i = [86 91 95 97 100 110 115 120 150 500 1000];
step = .01;
if alt < Z_i(1) || alt>Z_i(length(Z_i))
    n_i_array = [];
    return;
end
%   Gas coefficients
alpha_i = [0; 0; 0; 0; -.4; -.25];
a_i = [0; 6.986e20; 4.863e20; 4.487e20; 1.7e21; 3.305e21];
b_i = [0; .75; .75; .87; .691; .5];
Q_i = [0; -5.809644e-4; 1.366212e-4; 9.434079e-5; -2.457369e-4];
q_i = [0; -3.416248e-3; 0; 0; 0];
U_i = [0; 56.90311; 86; 86; 86];
u_i = [0; 97; 0; 0; 0];
W_i = [0; 2.70624e-5; 8.333333e-5; 8.333333e-5; 6.666667e-4];
w_i = [0; 5.008765e-4; 0; 0; 0];
%   Gas Data
R = 8.31432e3;
phi = 7.2e11;
T_7 = 186.8673;
T_11 = 999.2356;
%   Molecular Weight & Number Density based on values at 86 km & 500 km for
%   Hydrogen
n_i_86 = [1.129794e20; 8.6e16; 3.030898e19; 1.3514e18; 7.5817e14; 8e10];
n_i_alt = n_i_86;
sum_n = [ones(3,1)*n_i_86(1);ones(2,1)*sum(n_i_86(1:3));sum(n_i_86(1:5))];
M_i = [28.0134; 15.9994; 31.9988; 39.948; 4.0026; 1.00797];
M_0 = 28.9644;
n_int = zeros(size(n_i_86));
j = 1;
n_i_array = zeros(floor((alt-86)/division)+1,6);
for i = 1 : length(Z_i)-1
    if alt > Z_i(i)
        Z_start = Z_i(i);
        if alt > Z_i(i+1)
            Z_end = Z_i(i+1);
        else
            Z_end = alt;
        end
        for Z_0 = Z_start:step:Z_end-step
            Z_1 = Z_0+step;
            if Z_1 <= Z_i(5)
                M = ones(size(M_i))*M_0;
            else
                M = [(n_i_alt(1)*M_i(1))./sum_n(1:3);...
                    sum((n_i_alt(1:3).*M_i(1:3)))./sum_n(4:5);...
                    sum((n_i_alt(1:5).*M_i(1:5)))./sum_n(6)];
            end
            sum_n = [ones(3,1)*n_i_alt(1);ones(2,1)*sum(n_i_alt(1:3));sum(n_i_alt(1:5))];
            n_int = f_n(a_i,alpha_i,b_i,M,M_i,n_int,phi,...
                Q_i,q_i,R,sum_n,U_i,u_i,W_i,w_i,Z_i,Z_0,Z_1);
            n_i_alt(1:5) = n_i_86(1:5)*T_7/atmo_temp(Z_1).*exp(-n_int(1:5));
            if Z_1 < Z_i(9)
                n_i_alt(6) = 0;
            else
                tau = int_tau(alt);
                n_i_alt(6) = (T_11/atmo_temp(Z_1))^(1+alpha_i(6))*...
                    (n_i_86(6)*exp(-tau)-n_int(6));
            end
            if mod(Z_0,division) == 0
                n_i_array(j,:) = n_i_alt';
                j = j+1;
            end
                
        end
    end
end
n_i_end(1:5) = n_i_86(1:5)*T_7/atmo_temp(alt).*exp(-n_int(1:5));
if alt < Z_i(9)
    n_i_end(6) = 0;
else
    tau = int_tau(alt);
    n_i_end(6) = (T_11/atmo_temp(Z_1))^(1+alpha_i(6))*...
        (n_i_86(6)*exp(-tau)-n_int(6));
end
n_i_array(j,:) = n_i_end;
end

%% Pressure Function

function P = atmo_p(alt, T, sum_n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Program:    Atmospheric Pressure Calculation
%   Author:     Brent Lewis(RocketLion@gmail.com)
%               University of Colorado-Boulder
%   History:    Original-1/10/2007
%               Revision-1/12/2007-Corrected for changes in Matlab versions
%               for backward compatability-Many thanks to Rich
%               Rieber(rrieber@gmail.com)
%   Input:      alt:    Geometric altitude vector of desired pressure data[km]
%               T:      Temperature vector at given altitude points
%                       Required only for altitudes greater than 86 km[K]
%               sum_n:  Total number density of atmospheric gases[1/m^3]
%   Output:     P:      Pressure vector[Pa]
%   Note:       Must compute altitudes below 86 km and above 86 km on two
%               different runs to allow line up of altitudes and
%               temperatures
%   Examples:   atmo_p(0) = 101325 Pa
%               atmo_p(0:10) = Pressures between 0-10 km at 1 km increments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1
    T = [];
    sum_n = [];
end
%   Constants
N_A = 6.022169e26;
g_0 = 9.80665;
M_0 = 28.9644;
R = 8.31432e3;
r_E = 6.356766e3;
%   Geopotential/Geometric Altitudes used for Geometric Altitudes < 86 km
H = [0 11 20 32 47 51 71 84.852];
Z = r_E*H./(r_E-H);
Z(8) = 86;
%   Defined temperatures/lapse rates/pressures/density at each layer
T_M_B = [288.15 216.65 216.65 228.65 270.65 270.65 214.65];
L = [-6.5 0 1 2.8 0 -2.8 -2]/1e3;
P_ref = [1.01325e5 2.2632e4 5.4748e3 8.6801e2 1.1090e2 6.6938e1 3.9564];
%   Preallocation of Memory
P = zeros(size(alt));
for i = 1 : length(alt)
    Z_i = alt(i);
    if isempty(sum_n)
        index = find(Z>=Z_i)-1+double(Z_i==0);
        index = index(1);
        Z_H = r_E*Z_i/(r_E+Z_i);
        if L(index) == 0
            P(i) = P_ref(index)*exp(-g_0*M_0*(Z_H-H(index))*1e3/(R*T_M_B(index)));
        else
            P(i) = P_ref(index)*(T_M_B(index)/...
                (T_M_B(index)+L(index)*(Z_H-H(index))*1e3))^...
                (g_0*M_0/(R*L(index)));
        end
    else
        P(i) = sum_n(i)*R*T(i)/N_A;
    end
end
end

%% Temp Function

function Temp = atmo_temp(alt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Program:    Atmospheric Temperature Calculation
%   Author:     Brent Lewis(RocketLion@gmail.com)
%               University of Colorado-Boulder
%   History:    Original-1/09/2007
%   Input:      alt:    Geometric Altitude of desired altitude[scalar][km] 
%   Output:     Temp:   Temperature at desired altitude using values from 
%                       1976 Standard Atmosphere[K]
%   Note:       Only Valid Between 0 km and 1000 km
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Constants
r_E = 6.356766e3;
epsilon = 1e5*eps;
%   Defined temperatures at each layer
T = [288.15 216.65 216.65 228.65 270.65 270.65 ...
    214.65 186.95 186.8673 240 360 1000];
L = [-6.5 0 1 2.8 0 -2.8 -2 0 0 12 0];
%   Geopotential/Geometric Altitudes used for Geometric Altitudes < 86 km
H = [0 11 20 32 47 51 71];
Z = r_E*H./(r_E-H);
%   Geometric Altitudes used for Altitudes >86 km
Z(8:12) = [86 91 110 120 1000];
if alt < Z(1) || alt > (Z(12)+epsilon)
    error('Altitude must be 0-1000 km')
end
%   Temperature Calculation with Molecular Temperature below 86 km and
%   Kinetic Temperature above
if alt >= Z(1) && alt <= Z(8)
    Temp = interp1(Z,T,alt);
elseif alt > Z(8) && alt <= Z(9)
    Temp = T(9);
elseif alt > Z(9) && alt <= Z(10)
    a = 19.9429;
    A = -76.3232;
    T_c = 263.1905;
    Temp = T_c+A*sqrt(1-((alt-Z(9))/a)^2);
elseif  alt > Z(10) && alt <= Z(11)
    Temp = interp1(Z,T,alt);
elseif alt > Z(11)
    lambda = L(10)/(T(12)-T(11));
    xi = (alt-Z(11))*(r_E+Z(11))/(r_E+alt);
    Temp = T(12)-(T(12)-T(11))*exp(-lambda*xi);
end
end

%% Function 
function n_int = f_n(a_i,alpha_i,b_i,M,M_i,n_int,phi,...
    Q_i,q_i,R,sum_n,U_i,u_i,W_i,w_i,Z_i,Z_0,Z_1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Program:    Gas Integral Program
%   Author:     Brent Lewis(RocketLion@gmail.com)
%               University of Colorado-Boulder
%   History:    Original-1/10/2007
%   Input:      As defined in 1976 Standard Atmosphere
%   Output:     n_int:  Integral values computed using 5-point Simpsons
%                       Rule
%   Note:       Created for running in Atmospheric Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Constants
g_0 = 9.80665;
r_E = 6.356766e3;
T_c = 263.1905;
A_8 = -76.3232;
a_8 = 19.9429;
L_K_9 = 12;
T_7 = 186.8673;
T_9 = 240;
T_10 = 360;
T_11 = 999.2356;
T_inf = 1000;
%   Molecular Diffusion Coeffiecients
K_7 = 1.2e2;
alt_j = linspace(Z_0,Z_1,5);
n_i = zeros(6,length(alt_j));
for j = 1 : length(alt_j)
    Z = alt_j(j);
    %     Temperature Values
    if Z < Z_i(2)
        T = T_7;
        dT_dZ = 0;
    elseif Z < Z_i(6)
        T = T_c+A_8*sqrt(1-((Z-Z_i(2))/a_8)^2);
        dT_dZ = -A_8/a_8^2*(Z-Z_i(2))*(1-((Z-Z_i(2))/a_8)^2)^-.5;
    elseif Z < Z_i(8)
        T = T_9+L_K_9*(Z-Z_i(6));
        dT_dZ = L_K_9;
    elseif Z >= Z_i(8)
        lambda = L_K_9/(T_inf-T_10);
        xi = (Z-Z_i(8))*(r_E+Z_i(8))/(r_E+Z);
        T = T_inf-(T_inf-T_10)*exp(-lambda*xi);
        dT_dZ = lambda*(T_inf-T_10)*((r_E+Z_i(8))/(r_E+Z))^2*exp(-lambda*xi);
    end
    %     K Values
    if Z < Z_i(3)
        K = K_7;
    elseif Z < Z_i(7)
        K = K_7*exp(1-400/(400-(Z-95)^2));
    elseif Z >= Z_i(7)
        K = 0;
    end
    %     Gravity
    g = g_0*(r_E/(r_E+Z))^2;
    %     N
    if Z <= Z_i(5)
        M_N2 = M(1);
    else
        M_N2 = M_i(1);
    end
    n_i(1,j) = M_N2*g/(T*R)*1e3;
    %     O O2 Ar He
    D = a_i(2:5).*exp(b_i(2:5).*log(T/273.15))./sum_n(2:5);
    if K ~= 0
        f_Z = (g/(R*T)*D./(D+K).*(M_i(2:5)+M(2:5)*K./D+...
            alpha_i(2:5)*R/g*dT_dZ/1e3))*1e3;
    else
        f_Z = (g/(R*T)*(M_i(2:5)+alpha_i(2:5)*R/g*dT_dZ/1e3))*1e3;
    end
    if Z <= Z_i(4)
        vdk = Q_i(2:5).*([Z;Z;Z;Z]-U_i(2:5)).^2.*exp(-W_i(2:5).*...
            ([Z;Z;Z;Z]-U_i(2:5)).^3)+q_i(2:5).*(u_i(2:5)-[Z;Z;Z;Z]).^2.*...
            exp(-w_i(2:5).*(u_i(2:5)-[Z;Z;Z;Z]).^3);
    else
        vdk = Q_i(2:5).*([Z;Z;Z;Z]-U_i(2:5)).^2.*exp(-W_i(2:5).*...
            ([Z;Z;Z;Z]-U_i(2:5)).^3);
    end
    n_i(2:5,j) = f_Z+vdk;
    %     H
    if Z_1 < 150 || Z_1 > 500
        n_i(6,j) = 0;
    else
        D_H = a_i(6)*exp(b_i(6)*log(T/273.15))/sum_n(6);
        n_i(6,j) = phi/D_H*(T/T_11)^(1+alpha_i(6));
    end
end
h = alt_j(2)-alt_j(1);
n_int = n_int+(2*h/45*(7*n_i(:,1)+32*n_i(:,2)+12*n_i(:,3)+32*n_i(:,4)+7*n_i(:,5)));
end

%% Function
function TAU = int_tau(Z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Program:    Tau Integral Computation for Hydrogen Composition Program
%               in Atmospheric Model
%   Author:     Brent Lewis(RocketLion@gmail.com)
%               University of Colorado-Boulder
%   History:    Original-1/10/2007
%   Input:      Z:      Altitude value
%   Output:     TAU:    Integral Value
%   Note:       This program computes the value of Tau directly with the
%               integral done by hand and only the second integration limit
%               needing to be inputed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Constants
L_K_9 = 12;
T_10 = 360;
T_inf = 1000;
Z_10 = 120;
g_0 = 9.80665;
r_E = 6.356766e3;
R = 8.31432e3;
lambda = L_K_9/(T_inf-T_10);
M_H = 1.00797;
%   Value of Integration limit computed previously
tau_11 = 8.329503912749350e-004;
tau_Z = M_H*g_0*r_E^2/R*...
    log((exp(lambda*(Z-Z_10)*(r_E+Z_10)/(r_E+Z))-1)*T_inf+T_10)/...
    (lambda*T_inf*(r_E+Z_10)^2);
TAU = tau_Z-tau_11;
end

