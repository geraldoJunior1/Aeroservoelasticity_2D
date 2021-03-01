%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     SCRIPT PARA A DETERMINAÇÃO DE RESPOSTA      %
%               AEROSERVOELASTICA                 %
%                versao alpha-0.3                 %
%                                                 %
%       Desenvolvido por: Gerlado Junior          %
%               Harpia Aerodesign                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all 
close all

vel = 250; %velocidade
s = 10.25; % semi span
c = 2.5; % chord
a1 = 2*pi; %coef angular
rho = 1.225; % Densidade
muu = 250; % unidade de massa / area da asa
kappa_freq = 4; %flapping freq
theta_freq = 12.5; %pitch freq
xcm = 0.5*c; %posicao do CG
xf = 0.35*c; % posição do eixo elastico
dot_M_theta = -1.2; % correcao do termo de amortecimento não estacionario
e = xf/c - 0.25; % ecentricidade (CE - CA)
damping_Y_N = 1; %amortecimento

E_E = 0.15; % porcentagem da corda que a sup de controle ocupa
range = 5; %range de analise da sup
range = range * pi / 180;


seu_madruga = 0.25;
chirp_in = 1; % chirp inicial
chrip_out = 20; % chirp final

vento = rho*vel*c*s*[s/4 c/2]';

M11 = (muu*s^3*c)/3 ; 
M22 = muu*s*(c^3/3 - c*c*xf + xf*xf*c);
M12 = muu*s*s/2*(c*c/2 - c*xf );
M21 = M12;
M = [M11,M12; M21,M22];


k1 = (kappa_freq*pi*2)^2*M11; %kappa
k2 = (theta_freq*pi*2)^2*M22; %theta
E = [k1 0; 0 k2];

%matriz de rigidez
K = (rho*vel^2*[0,c*s^2*a1/4; 
                0,-c^2*s*e*a1/2])+[k1,0;0,k2] ;

%matriz de amortecimento

C = rho*vel*[c*s^3*a1/6,0; 
            -c^2*s^2*e*a1/4,-c^3*s*dot_M_theta/8];

%sups de controle
A_c = a1/pi*(acos(1-2*E_E) + 2*sqrt(E_E*(1-E_E))); 
B_c = -a1/pi*(1-E_E)*sqrt(E_E*(1-E_E));

forca_no_cont = rho*vel^2*c*s*[-s*A_c/4 c*B_c/2]';

%input de integração pro simulindo
M_K = inv(M)*K;
M_C = inv(M)*C;
M_FG = inv(M)*vento;
M_FC = inv(M)*forca_no_cont;


tmin = 0;  % start time
croc = 0.001; %passo
jacare = 10; % end time

tempoo = [0:croc:jacare]';

%sinal do sup de controle

t_fim = jacare * seu_madruga;
supp = zeros(size(tempoo)); 
passin_do_abencoado = sum(tempoo < t_fim);

for ii = 1:passin_do_abencoado
    supp(ii) = range*sin(2*pi*(chirp_in + (chrip_out -chirp_in)*ii/(2*passin_do_abencoado))*tempoo(ii));
end


%Rajada

rajada = zeros(size(tempoo));
rajada_dif = 0; % rajda "1 - cos" [muu/s]
tempo_de_rajada = 0.05;


lokao = jacare * tempo_de_rajada;
r_J = sum(tempoo < lokao); 
for ii = 1:r_J
    rajada(ii) = rajada_dif/2 * (1 - cos(2*pi*tempoo(ii)/lokao));
end


turbulencia_vertical_vel = 0;
tempode_de_turb = 1; 

t_fim = jacare * tempode_de_turb;
tempode_de_turb = sum(tempoo < t_fim); 
freq_max_da_turb = 20;
gravial = max(size(tempoo));

if rem(max(gravial),2) ~= 0
    gravial = gravial - 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


metade_grav = gravial / 2;
metade_grav_mais1 = gravial/2 + 1;
difff = 1/(gravial*croc);
n_freq = fix(freq_max_da_turb / difff ) + 1;  

for ii = 1:n_freq 
    a(ii) = 2 * rand - 1; % RE
    b(ii) = sqrt(1 - a(ii)*a(ii)) * (2*round(rand) - 1); %IM
end


c_f_r = (a + j*b); c_f_r(n_freq+1 : metade_grav_mais1) = 0;
c_f_r(metade_grav_mais1+1 : gravial) = conj(c_f_r(metade_grav:-1:2)); 
sup_turb = turbulencia_vertical_vel * real(ifft(c_f_r ));

for ii = 1:gravial
    rajada(ii) = rajada(ii) + sup_turb(ii);
end

% SIMULINDO
inrajada = [tempoo,rajada]; 
incontrole = [tempoo,supp]; 
[tout] = sim('malha');

dot_x_1 = saida1(:,1)*180/pi;
dot_x_2 = saida1(:,2)*180/pi;

x_1 = saida2(:,1)*180/pi; %kappa 
x_2 = saida2(:,2)*180/pi; %theta 

%plots
figure(1); 
plot(tempoo,supp,tempoo,rajada)
xlabel('t(s)'); 
ylabel('Angulo de deflexao(º) e Velocidade de rajada (m/s)')

figure(2) ; 
plot(tempoo,x_1,'r',tempoo,x_2,'b')
legend('flap', 'pitch');
xlabel('t(s)'); 
ylabel('angulo de controle (º)') 

figure(3); 
plot(tempoo,dot_x_1,'r',tempoo,dot_x_2,'b')
legend('flap', 'pitch');
xlabel('t(s)'); 
ylabel('Taxa de rotação (deg/s)')