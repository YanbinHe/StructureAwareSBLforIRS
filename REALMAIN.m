clc
clear

%% simluation on the SBL based IRS channel estimation with angle spread
% predefined values
ms_ante = 6; % the number of antennas at user equippment side
bs_ante = 16; % the number of antennas at bs
irs_ele = 16^2; % the number of irs elements
Res1 = 18; % the number of atoms in the dicitionay
timeslots = 1e5;
phase_offset = 0;

overheadp = 6; % the number of overhead w.r.t pilot signals
overheadi = [4,6,8,10]; % the number of overhead w.r.t reflecting pattern

angle_spread = 3; % the spread of angles in the domain of cos
angle = linspace(-1,1-2/Res1,Res1); 

SNRl = [5,10,15,20,25,30]; % in dB
SNRl_10 = 10.^(SNRl/10);

thres = 1e-4;
global thres_inner;
thres_inner = 1e-3;

numItr = 150; % the number of alternative optimization iterations
klevel = 50;

% codebook for optimal MIMO detection: search
symbolmatrix = pskmod([0:4-1], 4, phase_offset);
%% simulation settings
II = 4;
JJ = 6;
AVG = 100;% the number of trials

% results
errorl11 = zeros(II,JJ,AVG);
errorl1p = zeros(II,JJ,AVG);
errorl2 = zeros(II,JJ,AVG);
errorl3 = zeros(II,JJ,AVG);
errorc = zeros(II,JJ,AVG);
er = zeros(II,JJ,AVG);
r11 = zeros(II,JJ,AVG);
r12 = zeros(II,JJ,AVG);
r13 = zeros(II,JJ,AVG);
r2 = zeros(II,JJ,AVG);
s11 = zeros(II,JJ,AVG);
s12 = zeros(II,JJ,AVG);
s13 = zeros(II,JJ,AVG);
s2 = zeros(II,JJ,AVG);
sperf = zeros(II,JJ,AVG);
ro = zeros(II,JJ,AVG);
so = zeros(II,JJ,AVG);
time = zeros(3,II,JJ);
time1 = zeros(3,II,JJ);
time2 = zeros(3,II,JJ);
time3 = zeros(3,II,JJ);

X = 1/sqrt(ms_ante)*pilot_gen(ms_ante,max(overheadp));
irs_pattern = 1/sqrt(irs_ele)*((rand(irs_ele,max(overheadi)) > 0.5)*2-1); 
%%
for avg = 1 : AVG % for AVG trials
% different channel realization, stay constant within the coherence time
AoD_ms = randsample([1:Res1], 1);
g1 = zeros(Res1,1);
g1(AoD_ms) = 1;

Smpl = randsample([0:Res1-angle_spread], 1);
AoA_irs = mod(Smpl:1:Smpl+angle_spread-1,Res1)+1;
ga = zeros(Res1,1);
ga(AoA_irs) = 1;

AoD_irs = randsample([1:Res1], 1);
gd = zeros(Res1,1);
gd(AoD_irs) = 1;

Smpl = randsample([0:Res1-angle_spread], 1);
AoA_bs = mod(Smpl:1:Smpl+angle_spread-1,Res1)+1;  
g2 = zeros(Res1,1);
g2(AoA_bs) = 1;

% generate the true support, for the support recovery computation
gi = zeros(Res1,1);
indexi = [(AoA_irs(1)-AoD_irs(end)):1:(AoA_irs(end)-AoD_irs(1))] + Res1;
indexmodi = find(indexi > Res1);
indexi(indexmodi) = indexi(indexmodi) - Res1;
gi(indexi) = 1; 
gi = (fliplr(gi'))';

suppTrue = kron(gi,kron(g1,g2));
%%
% dictionary generation
A2 = generate_dict(bs_ante,Res1);
A_irs_d = generate_dict(irs_ele,Res1);
A_irs_a = generate_dict(irs_ele,Res1);
A1 = generate_dict(ms_ante,Res1); 

beta_1 = 1; % from ms to irs
beta_2 = 1; % from irs to bs
% path gain CN(0,1)
alpha_aco = (1*randn(angle_spread,1) + 1i*randn(angle_spread,1))/sqrt(2); % from ms to irs
alpha_a = zeros(Res1,1);
alpha_a(AoA_irs) = alpha_aco;

alpha_2co = (1*randn(angle_spread,1) + 1i*randn(angle_spread,1))/sqrt(2); % from irs to bs
alpha_2 = zeros(Res1,1);
alpha_2(AoA_bs) = alpha_2co;

H2 = sqrt(beta_2)*sqrt(bs_ante*irs_ele/angle_spread)*A2*(g2.*alpha_2)*gd'*A_irs_d';
H1 = sqrt(beta_1)*sqrt(ms_ante*irs_ele/angle_spread)*A_irs_a*(ga.*alpha_a)*g1'*A1';

% make sure that somehow the irs should point to the receiver, otherwise
% the received signal is too weak. But this is only for the simulation.
% This is not the case when it comes to real world.
for i = 1:max(overheadi)
    Htrue(:,:,i) = vec(H2*diag(irs_pattern(:,i))*H1);
    while norm(Htrue(:,:,i),'fro')<1e-5
        irs_pattern(:,i) =  1/sqrt(irs_ele)*(rand(irs_ele,1) > 0.5)*2-1; 
        Htrue(:,:,i) = vec(H2*diag(irs_pattern(:,i))*H1);
    end
end
%%
% generate signal
y_bar = kr(X.'*H1.',H2)*irs_pattern;
signal_power = norm(vec(y_bar))^2/length(vec(y_bar)); %average signal power per asymbol

for jj = 1:JJ % snr

SNR_10 = SNRl_10(jj);
noise_var = (signal_power)/SNR_10;

for ii = 1:II % overhead
% SNR and overhead settings
K1 = overheadp;
K2 = overheadi(ii);

x_pilot = X(:,1:K1);
IRS = irs_pattern(:,1:K2);
% different irs reflection patterns, the number is equal to #overheads,
H_p1 = IRS.'*kr(A_irs_a.',A_irs_d').';
H_p1 = H_p1(:,1:Res1);
H_p1_ori = H_p1;
H_p2 = x_pilot.'*conj(A1);
H_p2_ori = H_p2;
H = kron(kron(H_p1,H_p2),A2);
%% SNR part
% received noisy signals
y_tilde = vec(y_bar(1:K1*bs_ante,1:K2));
% generate noise
noise = sqrt(noise_var / 2)*(randn(size(y_tilde))+1i*randn(size(y_tilde)));
y = y_tilde+noise;
%% part 2: channel estimation with different techniques and compute the symbol error rate
% true channel for each IRS pattern
for i = 1:K2
    Htrue(:,:,i) = vec(H2*diag(IRS(:,i))*H1);
end

% transmitted data symbol
% generate bits stream
inPut = randi([0, 4-1], 1, timeslots);
Xun = pskmod(inPut, 4, phase_offset);
%% different techniques
%% SBL
[l2,timec,~,gc,Hrec] = classicSBL(numItr,H,Res1,noise_var,y,K2,IRS,A_irs_a,A_irs_d,A1,A2,H1,H2,SNRl(jj),thres);
% metrics compute
errorc(ii,jj,avg) = l2;
time(2,ii,jj) = time(2,ii,jj) + timec;

r2(ii,jj,avg) = recover_rate(suppTrue,gc);
se2 = ser_compute(Hrec,Htrue,SNRl(jj),Xun,symbolmatrix,ms_ante,bs_ante,K2);
s2(ii,jj,avg) = se2;
%% SVD
[l111,time111,~,gl1,Hrek1] = kroSBL1(noise_var,Res1,numItr,H,H_p1,H_p2,A1,A2,y,K2,IRS,A_irs_a,A_irs_d,H1,H2,SNRl(jj),thres,K2);
% metrics compute
errorl11(ii,jj,avg) = l111;
time1(1,ii,jj) = time1(1,ii,jj) + time111;

r11(ii,jj,avg) = recover_rate(suppTrue,gl1);
se11 = ser_compute(Hrek1,Htrue,SNRl(jj),Xun,symbolmatrix,ms_ante,bs_ante,K2);
s11(ii,jj,avg) = se11;
%% Alternating
[l12,time12,~,gl2,Hrek2] = kroSBL2(noise_var,Res1,numItr,H,H_p1,H_p2,A1,A2,y,K2,IRS,A_irs_a,A_irs_d,H1,H2,SNRl(jj),thres,K2);
errorl2(ii,jj,avg) = l12;
time2(1,ii,jj) = time2(1,ii,jj) + time12;

r12(ii,jj,avg) = recover_rate(suppTrue,gl2);
se12 = ser_compute(Hrek2,Htrue,SNRl(jj),Xun,symbolmatrix,ms_ante,bs_ante,K2);
s12(ii,jj,avg) = se12;
%% suboptimal
[l13,time13,~,gl3,Hrek3] = kroSBL3(noise_var,Res1,numItr,H,H_p1,H_p2,A1,A2,y,K2,IRS,A_irs_a,A_irs_d,H1,H2,SNRl(jj),thres,K2);
errorl3(ii,jj,avg) = l13;
time3(1,ii,jj) = time3(1,ii,jj) + time13;

r13(ii,jj,avg) = recover_rate(suppTrue,gl3);
se13 = ser_compute(Hrek3,Htrue,SNRl(jj),Xun,symbolmatrix,ms_ante,bs_ante,K2);
s13(ii,jj,avg) = se13;
%% perfect CSI 
sperf(ii,jj,avg) = ser_compute(Htrue,Htrue,SNRl(jj),Xun,symbolmatrix,ms_ante,bs_ante,K2);
%% OMP
[erroro,timeo,go,Hreo] = OMP(Res1,H,klevel,IRS,A_irs_a,A_irs_d,A1,A2,H1,H2,y,K2);
er(ii,jj,avg) = erroro;
time(3,ii,jj) = time(3,ii,jj) + timeo;
ro(ii,jj,avg) = recover_rate(suppTrue,go);
so(ii,jj,avg) = ser_compute(Hreo,Htrue,SNRl(jj),Xun,symbolmatrix,ms_ante,bs_ante,K2);
end
end
end
%%
save result_generate_by_user.mat
%% Plot part
%% ii-overhead jj-snr
results_kro1 = zeros(II,JJ);
results_kro2 = zeros(II,JJ);
results_kro3 = zeros(II,JJ);
results_cla = zeros(II,JJ);
results_omp = zeros(II,JJ);
results_r11 = zeros(II,JJ);
results_r12 = zeros(II,JJ);
results_r13 = zeros(II,JJ);
results_r2 = zeros(II,JJ);
results_ro = zeros(II,JJ);
results_s11 = zeros(II,JJ);
results_s12 = zeros(II,JJ);
results_s13 = zeros(II,JJ);
results_s2 = zeros(II,JJ);
results_sp = zeros(II,JJ);
results_so = zeros(II,JJ);

for ii = 1:II % for each overhead situation
    for jj = 1:JJ % snr
        % error
        results_kro1(ii,jj) = sum(errorl1(ii,jj,:))/AVG;
        results_kro2(ii,jj) = sum(errorl2(ii,jj,:))/AVG;
        results_kro3(ii,jj) = sum(errorl3(ii,jj,:))/AVG;
        results_cla(ii,jj) = sum(errorc(ii,jj,:))/AVG;
        results_omp(ii,jj) = sum(er(ii,jj,:))/AVG;
        % recovery rate
        results_r11(ii,jj) = sum(r11(ii,jj,:))/AVG;
        results_r12(ii,jj) = sum(r12(ii,jj,:))/AVG;
        results_r13(ii,jj) = sum(r13(ii,jj,:))/AVG;
        results_r2(ii,jj) = sum(r2(ii,jj,:))/AVG;
        results_ro(ii,jj) = sum(ro(ii,jj,:))/AVG;
        % ser
        results_s11(ii,jj) = sum(s11(ii,jj,:))/AVG;
        results_s12(ii,jj) = sum(s12(ii,jj,:))/AVG;
        results_s13(ii,jj) = sum(s13(ii,jj,:))/AVG;
        results_s2(ii,jj) = sum(s2(ii,jj,:))/AVG;
        results_sp(ii,jj) = sum(sperf(ii,jj,:))/AVG;
        results_so(ii,jj) = sum(so(ii,jj,:))/AVG;
    end
end
%%
figure('Units','centimeter','Position',[5 5 40 40]);
fontsizeman = 20;
ha = tight_subplot(2,2,[.1 .1],[.07 .07],[0.07 0.07]);
axes(ha(1)); 

snr1 = plot(SNRl,results_kro1(1,:),'-.o','Color',[0, 0.4470, 0.7410],'DisplayName','KroSBL-SVD');
hold on
snr2 = plot(SNRl,results_kro2(1,:),'-.x','Color',[0.6350, 0.0780, 0.1840],'DisplayName','KroSBL-Alternating');
hold on
snr3 = plot(SNRl,results_kro3(1,:),'-.m*','DisplayName','KroSBL-Appro');
hold on
snr4 = plot(SNRl,results_cla(1,:),'-.^','Color',[0.9290, 0.6940, 0.1250],'DisplayName','SBL');
hold on
snr5 = plot(SNRl,results_omp(1,:),'-.v','Color',[0.4940, 0.1840, 0.5560],'DisplayName','OMP');
hold on

plot(SNRl,results_kro1(4,:),'-o','Color',[0, 0.4470, 0.7410]);
hold on
plot(SNRl,results_kro2(4,:),'-x','Color',[0.6350, 0.0780, 0.1840]);
hold on
plot(SNRl,results_kro3(4,:),'-m*');
hold on
plot(SNRl,results_cla(4,:),'-^','Color',[0.9290, 0.6940, 0.1250]);
hold on
plot(SNRl,results_omp(4,:),'-v','Color',[0.4940, 0.1840, 0.5560]);
hold on

% legend([snr1,snr2,snr3,snr4,snr5],'SVD-KroSBL', 'AM-KroSBL','Classic KroSBL [9,10]','Classic SBL','OMP');
% legend('boxoff')
title('(a)')
xlabel('SNR')
ylabel('NMSE')
grid on
set(gca, 'yscale', 'log');
set(gca, 'yscale', 'log');
set(0,'DefaultLineLineWidth',3)
set(0,'DefaultAxesFontSize',fontsizeman)
set(0,'DefaultLineMarkerSize',14)
set(0,'DefaultAxesFontWeight','bold')
set(gca,'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontSize',fontsizeman)
set(get(gca,'Ylabel'),'FontSize',fontsizeman)
set(get(gca,'Title'),'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontWeight','bold')
set(get(gca,'Ylabel'),'FontWeight','bold')
set(get(gca,'Title'),'FontWeight','bold')
%
% figure
% subplot(2,2,2)
axes(ha(2)); 
for ii = 1 % for each SNR case
srr1 = plot(SNRl,results_r11(ii,:),'-.o','Color',[0, 0.4470, 0.7410],'DisplayName','KroSBL-SVD');
hold on
srr2 = plot(SNRl,results_r12(ii,:),'-.x','Color',[0.6350, 0.0780, 0.1840],'DisplayName','KroSBL-Alternating');
hold on
srr3 = plot(SNRl,results_r13(ii,:),'-.m*','DisplayName','[12]');
hold on
srr4 = plot(SNRl,results_r2(ii,:),'-.^','Color',[0.9290, 0.6940, 0.1250],'DisplayName','SBL');
hold on
srr5 = plot(SNRl,results_ro(ii,:),'-.v','Color',[0.4940, 0.1840, 0.5560],'DisplayName','OMP');
hold on
end
for ii = 4 % for each SNR case
plot(SNRl,results_r11(ii,:),'-o','Color',[0, 0.4470, 0.7410],'DisplayName','KroSBL-SVD');
hold on
plot(SNRl,results_r12(ii,:),'-x','Color',[0.6350, 0.0780, 0.1840],'DisplayName','KroSBL-Alternating');
hold on
plot(SNRl,results_r13(ii,:),'-m*','DisplayName','[12]');
hold on
plot(SNRl,results_r2(ii,:),'-^','Color',[0.9290, 0.6940, 0.1250],'DisplayName','SBL');
hold on
plot(SNRl,results_ro(ii,:),'-v','Color',[0.4940, 0.1840, 0.5560],'DisplayName','OMP');
hold on
end
xlabel('SNR')
ylabel('SRR')
title('(b)')
grid on
% legend([srr1,srr2,srr3,srr4,srr5],'SVD-KroSBL', 'AM-KroSBL','Classic KroSBL [9,10]','Classic SBL','OMP');
% legend('boxoff')
% set(gca, 'yscale', 'log');
% set(gca, 'yscale', 'log');
set(0,'DefaultLineLineWidth',3)
set(0,'DefaultAxesFontSize',fontsizeman)
set(0,'DefaultLineMarkerSize',14)
set(0,'DefaultAxesFontWeight','bold')
set(gca,'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontSize',fontsizeman)
set(get(gca,'Ylabel'),'FontSize',fontsizeman)
set(get(gca,'Title'),'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontWeight','bold')
set(get(gca,'Ylabel'),'FontWeight','bold')
set(get(gca,'Title'),'FontWeight','bold')

%
% figure
% subplot(2,2,3)
axes(ha(3)); 
for ii = 1
plot(SNRl,reshape(time1(1,ii,:)/AVG,JJ,1),'-.o','Color',[0, 0.4470, 0.7410],'DisplayName','KroSBL with SVD');
hold on
plot(SNRl,reshape(time2(1,ii,:)/AVG,JJ,1),'-.x','Color',[0.6350, 0.0780, 0.1840],'DisplayName','KroSBL with Alternating');
hold on
plot(SNRl,reshape(time3(1,ii,:)/AVG,JJ,1),'-.m*','DisplayName','[12]');
hold on
plot(SNRl,reshape(time(2,ii,:)/AVG,JJ,1),'-.^','Color',[0.9290, 0.6940, 0.1250],'DisplayName','SBL');
hold on
plot(SNRl,reshape(time(3,ii,:)/AVG,JJ,1),'-.v','Color',[0.4940, 0.1840, 0.5560],'DisplayName','OMP');
hold on
end

dash = plot(SNRl,zeros(length(SNRl),1),'-.k');
solid = plot(SNRl,zeros(length(SNRl),1),'-k');

for ii = 4
t1 = plot(SNRl,reshape(time1(1,ii,:)/AVG,JJ,1),'-o','Color',[0, 0.4470, 0.7410],'DisplayName','KroSBL with SVD');
hold on
t2 = plot(SNRl,reshape(time2(1,ii,:)/AVG,JJ,1),'-x','Color',[0.6350, 0.0780, 0.1840],'DisplayName','KroSBL with Alternating');
hold on
t3 = plot(SNRl,reshape(time3(1,ii,:)/AVG,JJ,1),'-m*','DisplayName','[12]');
hold on
t4 = plot(SNRl,reshape(time(2,ii,:)/AVG,JJ,1),'-^','Color',[0.9290, 0.6940, 0.1250],'DisplayName','SBL');
hold on
t5 = plot(SNRl,reshape(time(3,ii,:)/AVG,JJ,1),'-v','Color',[0.4940, 0.1840, 0.5560],'DisplayName','OMP');
hold on
end
set(gca, 'yscale', 'log');
xlabel('SNR')
ylabel('Run Time (Seconds)')
grid on


legend([dash,solid],'$K_{\mathrm{I}}=4$','$K_{\mathrm{I}}=10$','Interpreter', 'LaTex');
legend('boxoff')

set(0,'DefaultLineLineWidth',3)
set(0,'DefaultAxesFontSize',fontsizeman)
set(0,'DefaultLineMarkerSize',14)
set(0,'DefaultAxesFontWeight','bold')
set(gca,'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontSize',fontsizeman)
set(get(gca,'Ylabel'),'FontSize',fontsizeman)
set(get(gca,'Title'),'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontWeight','bold')
set(get(gca,'Ylabel'),'FontWeight','bold')
set(get(gca,'Title'),'FontWeight','bold')
title('(c)')


%
% figure
% subplot(2,2,4)
axes(ha(4)); 
for ii = 1
ser1 = plot(SNRl,results_s11(ii,:),'-.','Color',[0, 0.4470, 0.7410],'DisplayName','KroSBL-SVD');
hold on
ser2 = plot(SNRl,results_s12(ii,:),'-.','Color',[0.6350, 0.0780, 0.1840],'DisplayName','KroSBL-Alternating');
hold on
ser3 = plot(SNRl,results_s13(ii,:),'-.m','DisplayName','[12]');
hold on
ser4 = plot(SNRl,results_s2(ii,:),'-.','Color',[0.9290, 0.6940, 0.1250],'DisplayName','SBL');
hold on
ser5 = plot(SNRl,results_so(ii,:),'-.','Color',[0.4940, 0.1840, 0.5560],'DisplayName','OMP');
hold on
end

for ii = 1
ser1 = plot(SNRl,results_s11(ii,:),'o','Color',[0, 0.4470, 0.7410],'DisplayName','KroSBL-SVD');
hold on
ser2 = plot(SNRl,results_s12(ii,:),'x','Color',[0.6350, 0.0780, 0.1840],'DisplayName','KroSBL-Alternating');
hold on
ser3 = plot(SNRl,results_s13(ii,:),'m*','DisplayName','[12]');
hold on
ser4 = plot(SNRl,results_s2(ii,:),'^','Color',[0.9290, 0.6940, 0.1250],'DisplayName','SBL');
hold on
ser5 = plot(SNRl,results_so(ii,:),'v','Color',[0.4940, 0.1840, 0.5560],'DisplayName','OMP');
hold on
end

for ii = 4
ser21 = plot(SNRl,results_s11(ii,:),'-o','Color',[0, 0.4470, 0.7410],'DisplayName','KroSBL-SVD');
hold on
ser22 = plot(SNRl,results_s12(ii,:),'-x','Color',[0.6350, 0.0780, 0.1840],'DisplayName','KroSBL-Alternating');
hold on
ser23 = plot(SNRl,results_s13(ii,:),'-m*','DisplayName','[12]');
hold on
ser24 = plot(SNRl,results_s2(ii,:),'-^','Color',[0.9290, 0.6940, 0.1250],'DisplayName','SBL');
hold on
ser25 = plot(SNRl,results_so(ii,:),'-V','Color',[0.4940, 0.1840, 0.5560],'DisplayName','OMP');
hold on
end


serp = plot(SNRl,results_sp(ii,:),':+','color',[0.4660, 0.6740, 0.1880],'DisplayName','Perfect CSI');
hold on
% serp = plot(SNRl,results_sp(ii,:),'+','color',[0.4660, 0.6740, 0.1880],'DisplayName','Perfect CSI');
% hold on

legend([ser1,ser2,ser3,ser4,ser5,serp],'SVD-KroSBL', 'AM-KroSBL','KroSBL [9,10]','Classic SBL','OMP','Perfect CSI','FontSize',fontsizeman,'FontWeight','normal');
legend('boxoff')

set(0,'DefaultLineLineWidth',3)
set(0,'DefaultAxesFontSize',fontsizeman)
set(0,'DefaultLineMarkerSize',14)
set(0,'DefaultAxesFontWeight','bold')
set(gca,'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontSize',fontsizeman)
set(get(gca,'Ylabel'),'FontSize',fontsizeman)
set(get(gca,'Title'),'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontWeight','bold')
set(get(gca,'Ylabel'),'FontWeight','bold')
set(get(gca,'Title'),'FontWeight','bold')
title('(d)')
xlabel('SNR')
ylabel('SER')
grid on
set(gca, 'yscale', 'log');