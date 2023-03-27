% =========================================================================
% Title       : Generate test plot
% File        : testplot.m
% -------------------------------------------------------------------------
% Revisions   :
%   Date       Version  Author  Description
%   11-dec-11  1.0      studer  initial version
% -------------------------------------------------------------------------
%   (C) 2006-2011 Communication Theory Group                      
%   ETH Zurich, 8092 Zurich, Switzerland                               
%   Author: Dr. Christoph Studer (e-mail: studer@rice.edu)     
% =========================================================================

% -- initialization
clear all
clc
maxfiles = 10;

% -- read and consolidate results
PCTC_3GPP_R13 = consolidate('ERR_PCTC_3GPP_V3_R13_640b_I10',maxfiles);
PCTC_3GPP_R12 = consolidate('ERR_PCTC_3GPP_V3_R12_640b_I10',maxfiles);


% -- BER
figure(1);
semilogy(PCTC_3GPP_R13.TxRx.Sim.EbNo_dB_list,PCTC_3GPP_R13.BER,'rs-')
hold on
semilogy(PCTC_3GPP_R12.TxRx.Sim.EbNo_dB_list,PCTC_3GPP_R12.BER,'kx-')
hold off
grid on
xlabel('E_b/N_o [dB]')
ylabel('BER')
axis([0 3 1e-6 1])
legend('3GPP R=1/3 L=640 I=10','3GPP R=1/2 L=640 I=10')

% -- FER
figure(2);
semilogy(PCTC_3GPP_R13.TxRx.Sim.EbNo_dB_list,PCTC_3GPP_R13.FER,'rs-')
hold on
semilogy(PCTC_3GPP_R12.TxRx.Sim.EbNo_dB_list,PCTC_3GPP_R12.FER,'kx-')
hold off
grid on
xlabel('E_b/N_o [dB]')
ylabel('FER')
axis([0 3 1e-4 1])
legend('3GPP R=1/3 L=640 I=10','3GPP R=1/2 L=640 I=10')

