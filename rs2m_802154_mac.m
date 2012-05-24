%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2012 Telecom SudParis.
%
% Permission is hereby granted, free of charge, to any person obtaining
% a copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to
% permit persons to whom the Software is furnished to do so, subject to
% the following conditions:
%
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
% SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB M-file rs2m_802154_mac.m
% Authors : Mohamed-Haykel Zayani and Vincent Gauthier
% Emails : {mohamed-haykel.zayani, vincent.gauthier}@it-sudparis.eu
% Address : Laboratory CNRS S.A.M.O.V.A.R. - Dept RS2M
% Telecom Sud Paris * 9 rue C. Fourier * 91011 EVRY CEDEX * FRANCE
% Created : April 28th, 2010
% Updated : May 18th, 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Brief Description:
% ------------------
% The m-file analyzes an Institute of Electrical and Electronics Engineers 
% (IEEE) 802.15.4 Machine Access Control (MAC) layer channel in which
% multiple nonsaturated stations compete for communicate with a sink. It is
% inspired from the IEEE 802.11 Model [2] developed by David Griffith and
% Michael Souryal (Emerging and Mobile Network Technologies Group,
% Information Technology Laboratory, National Institute of Standards and
% Technology). The objective behind this "adaptation" is to model a physical
% layer (PHY), including path loss and shadowing effects. The particularity
% of the proposed model lies in overstepping the node range disk shaped and
% taking into consideration the called “transitional area”. The model
% relies on the approach of Zuniga and Krishnamachari (references are in
% ZunPhyModel.m file). The function ZunPhyModel performs the calculations
% at the PHY level and determines the probability of good frame reception
% towards channel (signal-to-noise ratio) and radio (modulation and coding)
% setups.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% -----------
% [1] SGIP NIST Smart Grid Collaboration Site, PAP02: Wireless 
% Communicattions for the Smart Grid (6.1.5), Call for Input to Task6, 2010. 
% http://collaborate.nist.gov/twiki-sggrid/bin/view/SmartGrid/PAP02Wireless#Call_for_Input_to_Task_6. 
%
% [2] Park, P., Di Marco, P., Soldati, P., Fischione, C. and Johansson,
% K.H., "A Generalized Markov Chain Model for Effective Analysis of Slotted
% IEEE 802.15.4", Mobile Adhoc and Sensor Systems, IEEE 6th International
% Conference, Macau, pp. 130-139, 2009.
%
% [3] Pollin, S., Ergen, L., Ergen, S. C., Bougard, B., Catthoor, F.,
% Bahai,A., and Varaiya, P., “Performance analysis of slotted carrier sense
% IEEE 802.15.4 acknowledged uplink transmissions,” in Proc. of IEEE WCNC,
% pp. 1559–1564, 2008.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
% Get the time when we start computations:
start_time = clock;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define global parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of stations:
N_stations = 10;
% System size (buffer + service, so buffer holds K-1 frames):
K = 51;
% Set the MAC contention scheme:
access_schemes = {'basic'};
MAC_scheme = access_schemes{1};
%Physical layer data rate in bits/s, based on the IEEE 802.15.4 mode
data_rate = 19.2*10^3;
%The minimum value of the backoff exponent (BE) in the CSMA-CA algorithm
macMinBE=3;
%The maximum value of the backoff exponent (BE) in the CSMA-CA algorithm
macMaxBE=5;
%The maximum number of backoffs the CSMA-CA algorithm will attempt before
%declaring a channel access failure
macMaxCSMABackoffs=4;
%The maximum number of retries allowed after a transmission failure
macMaxFrameRetries=3;
% Initial retry backoff window size:
W0 = 8;
% Maximum number of frame transmission retries:
alpha = macMaxFrameRetries;
% Maximum number of backoffs, maximum backoff window size is Wmax = 2^m *
% W0:
m = macMaxBE-macMinBE;
% Compute the vector of backoff window sizes.  If the number of
% transmission retries is bigger than the maximum backoff window size, the
% window expands by a factor of 2 with each backoff until the limit m is
% reached and remains fixed at Wmax thereafter.
W(1) = W0;
if alpha > m,
    if m > 0,
        W(2:m+1) = 2.^(1:1:m) * W(1);
        W(m+1:alpha+1) = W(m+1) * ones(1, alpha - m + 1);
    elseif m == 0,
        W(2:alpha+1) = W(1) * ones(1,alpha);
    end
else
    W(2:alpha+1) = 2.^(1:1:alpha) * W(1);
end

% Size of MAC frame payload (Data Field), in bits: 
L_application = 121*8;
% Size of overhead added in PHY layer (Preamble + Start of Packet Delimiter
% + PHY Header), in bits:
L_overhead = 6*8;
% Size of frame payload in bits (application + overhead) MAX should be 127 bytes
% as defined in the standard
L_payload = L_application + L_overhead;
% ACK frame size in bits at PHY layer
L_ACK = 11*8;
% Mean propagation delay in seconds:
T_prop = 222e-9;

% Vector of lambda values where each entry is the per-station frame
% arrival rate in frames/s:
lvec = 0.5:0.1:7;
% Set the basic backoff time period used by the CSMA-CA algorithm, in
% seconds
aUnitBackoffPeriod = 16*20e-6;
% Idle state length without generating packets
L0 = 0;
% Coefficient used to translate the time length of a frame to slot length
% (bits/slot)
A = 80;

% IFS Definition in function of payload size  
if L_payload > 18*8
    t_IFS=40*16e-6;
else
    t_IFS=12*16e-6;
end
% The maximum time to wait for an acknowledgment frame to
% arrive following a transmitted data frame
macACKWaitDuration = 120*16e-6;
% RX-to-TX or TX-to-RX maximum turnaround time
aTurnaroudTime = 12*16e-6;
% The time used to detect if the channel is busy or not
sensing_time = 8*16e-6;

% C_t is the length of a collision slot or a slot in which the frame is
% lost due to bit errors, in seconds
C_t = (L_payload/data_rate) + T_prop + macACKWaitDuration;
% S_t is the length of a successful transmission cycle, in seconds.
S_t = (L_payload/data_rate) + aTurnaroudTime + aUnitBackoffPeriod + (L_ACK/data_rate) + 2*T_prop + t_IFS;

% Allocate memory for output arrays:
% Mean service time:
ET = zeros(size(lvec));
% Standard deviation of the service time:
Std_dev = zeros(size(lvec));
% Blocking probability:
P_blocking = zeros(size(lvec));
% Probability station is idle:
P_idle = zeros(size(lvec));
% Frame transmission failure probability:
P_failure = zeros(size(lvec));
% Access channel failure probability:
Pcf = zeros(size(lvec));
% Frame transmissions failure probability:
Pcr = zeros(size(lvec));
% Average number in system:
L_value = zeros(size(lvec));
% Reliability metric (Probability that application-generated packet is
% successfully sent):
Reliability = zeros(size(lvec));
% Average wait time (using Little's Theorem):
D_value = zeros(size(lvec));
% Average throughput:
S_avg = zeros(size(lvec));
% Instantaneous throughput:
S_inst = zeros(size(lvec));
% Channel access failure probability from (19) in Park
Pcf = zeros(size(lvec));
% Packet discarded due to retry limits probability from (20) in Park
Pcr = zeros(size(lvec));
% Alpha and beta probabilities from (17) and (18) in Park
Alpha = zeros(size(lvec));
Beta = zeros(size(lvec));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute P_e, the Frame Error Rate (FER):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nakagami fading parameter (m=1 is Rayleigh fading, m=inf is no fading):
% Nakagami_m = inf; 
% Shadowing standard deviation, in dB:
sigma_s = 4;
% The FER is the probability that the link is up, based on the channel
% conditions. ZunPhyModel is a function that returns the probability that a
% packet is successfully received towards channel conditions. The function
% relies on the model proposed by Zuniga & Krishnamachari (references 
% are in the ZunPhyModel.m file):
P_e = 1-ZunPhyModel(sigma_s, L_payload, data_rate);

% Tolerances for convergence of p0 and tau:
p0_tolerance = 1e-10;

% Initialize counter to point to the first element of the vector lvec:
counter = 1;

for lambda = lvec,
   
     % Get the time when we start computations for this value of lambda:
    loop_start = clock;
    
    % Initialize the values for p0 and p so that the while look will
    % execute on the first pass:
    p0 = -1;
    p = zeros(1,K+1);
    
    while abs(p0 - p(1)) > p0_tolerance,
        
        % Update p0 with the value that we got on the previous iteration:
        p0 = p(1);
        
        % Solve the four non-linear equations (16), (17) for alpha1 and 
        % alpha2 and (18) in Park and inspired from Pollin [6]. The 
        % equations are written in their long form to highlight the four 
        % probabilities to determine : tau (probability that a node 
        % attempts a first carrier sensing as z(1)), alpha1 (probability of
        % finding channel busy during CCA1 due to data transmission as 
        % z(2)), alpha2 (probability if finding channel busy during CCA1 
        % due to ACK transmission as z(3)) and beta (probability 
        % probability of finding channel busy during CCA2 as z(4))
                
                
        f=@(z)([
                % Equation (16) in Park
                z(1)-(((1-(z(2)+z(3)+(1-z(2)-z(3))*...
                z(4))^(macMaxCSMABackoffs+1))/...
                (1-(z(2)+z(3)+(1-z(2)-z(3))*z(4))))*...
                ((1-(1 - (1 - P_e) * (1 - (1 - (1-(1-p0)*z(1))^(N_stations-1))))*((z(2)+z(3)+(1-z(2)-z(3))*z(4))^(macMaxCSMABackoffs+1))^...
                (macMaxFrameRetries+1))/...
                (1-(1 - (1 - P_e) * (1 - (1 - (1-(1-p0)*z(1))^(N_stations-1))))*(1-(z(2)+z(3)+(1-z(2)-z(3))*z(4))^(macMaxCSMABackoffs+1))))*...
                ( (W0/2)*(1+2*(z(2)+z(3)+...
                (1-z(2)-z(3))*z(4)))*...
                (1+(1 - (1 - P_e) * (1 - (1 - (1-(1-p0)*z(1))^(N_stations-1))))*(1-(z(2)+z(3)+...
                (1-z(2)-z(3))*z(4))^(macMaxCSMABackoffs+1)))+S_t*(1-(z(2)+z(3)+...
                (1-z(2)-z(3))*z(4))^2)*...
                (1+(1 - (1 - P_e) * (1 - (1 - (1-(1-p0)*z(1))^(N_stations-1))))*(1-(z(2)+z(3)+(1-z(2)-z(3))*z(4))^(macMaxCSMABackoffs+1)))+...
                ((L0*p0)/(1-p0))*...
                (((1 - (1 - P_e) * (1 - (1 - (1-(1-p0)*z(1))^(N_stations-1))))*(1-(z(2)+z(3)+(1-z(2)-z(3))*z(4))^2))^2*...
                (((1 - (1 - P_e) * (1 - (1 - (1-(1-p0)*z(1))^(N_stations-1))))*(1-(z(2)+z(3)+(1-z(2)-z(3))*z(4))^2))^(macMaxFrameRetries-1)+1)+1) )^(-1)),
                
                % Equation (17) for alpha1 in Park
                z(2)-(L_payload/A)*((1-(1-(1-p0)*z(1))^(N_stations-1))*(1-z(2)-z(3))*(1-z(4))),
                    
                % Equation (17) for alpha2 in Park
                z(3)-(L_ACK/A)*(N_stations*(1-p0)*z(1)*((1-(1-p0)*z(1))^(N_stations-1))*(1-(1-(1-p0)*z(1))^(N_stations-1))...
                *(1-z(2)-z(3))*(1-z(4))*(1/(1-(1-(1-p0)*z(1))^N_stations))),
                
                % Equation (18) in Park
                z(4)-(( 1-((1-(1-p0)*z(1))^(N_stations-1))+N_stations*(1-p0)*z(1)*(1-(1-p0)*z(1))^(N_stations-1) )...
                / ( 2-(1-(1-p0)*z(1))^(N_stations)+N_stations*(1-p0)*z(1)*(1-(1-p0)*z(1))^(N_stations-1) ))
                ]);
                     
        % Initial solution for the system
        param0=[0.3,0.8,0.05,0.5];
                
        % Options of the solving method 'fsolve'
        options=optimset('MaxFunEvals',100000,'MaxIter',10000,'Display','off');
                
        % Results after solving the system
        out=fsolve(f,param0,options);
                
        
        % Assignment of found results 
        tau = out(1);
        alpha1_CCA1 = out(2);
        alpha2_CCA1 = out(3);
        alpha_CCA1 = alpha1_CCA1 + alpha2_CCA1;
        beta_CCA2 = out(4);
        
        % Computation of the probability of a collision.
        % This is the probability that at least one other station
        % transmits during the desired time slot:
        P_col = 1 - (1-(1-p0)*tau)^(N_stations-1);
                
        % Computation of the equivalent probability of a failed
        % transmission attempt.  This is the probability that a
        % sent frame is received correctly and it does not collide
        % with any other frames:
        P_fail = 1 - (1 - P_e) * (1 - P_col);
        
        % Computation of x and y parameters like in Park (section III, 
        % page 3)
        x = alpha_CCA1+(1-alpha_CCA1)*beta_CCA2;
        y = P_fail*(1-x^(macMaxCSMABackoffs+1));
     
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % With tau in hand, get throughput and other statistics.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Computing the average backoff period (approximation) based on
        % (35) in Park
        E_T = 0;
        for i = 0:macMaxCSMABackoffs
        E_T = E_T + 2*sensing_time + ( ((max(alpha_CCA1,(1-alpha_CCA1)*beta_CCA2)).^(i)) ./ ...
            (sum(max(alpha_CCA1,(1-alpha_CCA1)*beta_CCA2).^(0:macMaxCSMABackoffs) ) ) .* ...
            ( sum( (((W0*2.^(0:i))-1)/2)*aUnitBackoffPeriod+2*sensing_time.*(0:i) ) ));
        end
        % Computing the average delay for a successfully received packet
        % (approximation) based on (34) in Park
        ET(counter) = sum( (1-P_fail.*(1-x^(macMaxCSMABackoffs+1))).*P_fail.^(0:macMaxFrameRetries)...
            .*(1-x^(macMaxCSMABackoffs+1)).^(0:macMaxFrameRetries)...
                    *(1/(1-(P_fail*(1-x^(macMaxCSMABackoffs+1)))^(macMaxFrameRetries+1))) ...
                    .* (S_t+(0:macMaxFrameRetries)*C_t+((0:macMaxFrameRetries)+1).*E_T)  );
        
        Std_dev(counter) = sqrt( sum( (1-P_fail.*(1-x^(macMaxCSMABackoffs+1))).*P_fail.^(0:macMaxFrameRetries)...
            .*(1-x^(macMaxCSMABackoffs+1)).^(0:macMaxFrameRetries)...
                    *(1/(1-(P_fail*(1-x^(macMaxCSMABackoffs+1)))^(macMaxFrameRetries+1))) ...
                    .* (S_t+(0:macMaxFrameRetries)*C_t+((0:macMaxFrameRetries)+1).*E_T).^2) - ET(counter)^2);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute M/M/1/K queue model state probabilities
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % The utilization is the ratio of the arrival and service rates,
        % and the service rate is the reciprocal of the mean service time.
        rho = lambda * ET(counter);
        
        % p_0 is the reciprocal of the sum of powers of rho, and p_i is
        % given by p_i = rho^i * p_0:
        p(1) = 1/sum(rho.^(0:1:K));
        p(2:K+1) = rho.^(1:1:K) * p(1);

    end % End of p0-updating while loop
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute metrics associated with the current lambda value in the
    % vector lvec:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Blocking probability
    P_blocking(counter) = p(K+1);
    % Probability that the station is idle (this is p_0, the probability
    % that the queue is in state 0):
    P_idle(counter) = p(1);
    % Frame transmission failure probability (due to collisions and channel
    % errors):
    P_failure(counter) = P_fail;
    % Average number of frames in the system at each station:
    L_value(counter) = (0:K) * p';
    % Channel access failure probability from (19) in Park
    Pcf(counter) = (1/(1-y))*x^(macMaxCSMABackoffs+1)*(1-y^(alpha+1));
    % Packet discarded due to retry limits probability from (20) in Park
    Pcr(counter) = y^(alpha+1);
    % Alpha and beta probabilities
    Alpha(counter) = alpha1_CCA1+alpha2_CCA1;
    Beta(counter) = beta_CCA2;
    % Reliability (probability that frame is not blocked or lost):
    Reliability(counter) = (1-P_blocking(counter)) * (1 - Pcf(counter)) * (1 - Pcr(counter));
    % Average wait time (using Little's Theorem), in seconds:
    D_value(counter) = L_value(counter) / (lambda * (1-P_blocking(counter)));
    % Throughput (in bits/s):
    S_avg(counter) = lambda * Reliability(counter) * L_application;
    % Instantaneous throughput (in bits/s):
    S_inst(counter) = L_payload / ET(counter);
    
    X_c(counter)=x;
    Y_c(counter)=y;
    
    % Generate the status text block and display it in the Command
    % Window:
    out_s = [sprintf('lambda = %f\np0 is now %f\n',lambda,p0), ...
        sprintf('Got P_fail, P_e, P_col and tau; they are %f, %f, %f and %f\n',P_fail,P_e, P_col, tau), ...
        sprintf('Got alpha and beta; they are %f and %f\n',alpha_CCA1, beta_CCA2), ...
        sprintf('Getting state probabilities.\n'), ...
        sprintf('p_B = %f, p_f = %f, R = %f\n',P_blocking(counter),P_failure(counter),Reliability(counter))];
    disp(out_s);
    
    % Increment the lvec counter:
    counter = counter + 1;
    
    fprintf('Time for lambda = %f is %f seconds\n\n',...
        lambda,etime(clock,loop_start))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute total elapsed time and break it out into hours, minutes, and
% seconds.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First get the current time:
end_time = clock;
% The difference in the start and end times is returned in units of sec by
% etime():
elapsed_time = etime(end_time,start_time);
% There are 3600 seconds in an hour:
elapsed_hours = floor(elapsed_time/3600);
% Divide the remaining elapsed time by 60 and round down to get minutes:
elapsed_minutes = floor((elapsed_time - 3600*elapsed_hours)/60);
% The remainder is measured in seconds:
elapsed_seconds = elapsed_time - 3600*elapsed_hours - 60*elapsed_minutes;
fprintf('Total execution time is %d hours, %d minutes, %3.1f seconds\n',...
    elapsed_hours,elapsed_minutes,elapsed_seconds)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots of metrics evolution associated to lambda increment  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cnt = 1:counter-1;
figure,
hold on,
grid on,
plot(lvec*L_application,D_value(cnt),'-or');
title( ['Average Wait Time VS Offered Load for ',num2str(N_stations),' Nodes']);
xlabel('Offered Load (bits/application/node)');
ylabel('Average wait time (seconds)');
figure,
hold on,
grid on,
plot(lvec*L_application,Reliability(cnt),'-ob');
title( ['Reliability VS Offered Load for ',num2str(N_stations),' Nodes']);
xlabel('Offered Load (bits/application/node)');
ylabel('Reliability');
figure,
hold on,
grid on,
plot(lvec*L_application,P_failure(cnt),'-oy');
title( ['Failure Probability VS Offered Load for ',num2str(N_stations),' Nodes']);
xlabel('Offered Load (bits/application/node)');
ylabel('Failure Probability');
figure,
hold on,
grid on,
plot(lvec*L_application,Alpha(cnt),'-*');
hold on,
plot(lvec*L_application,Beta(cnt),'-o');
title( ['Alpha & Beta Probabilities VS Offered Load for ',num2str(N_stations),' Nodes']);
legend('Alpha','Beta');
xlabel('Offered Load (bits/application/node)');
ylabel('Alpha/Beta');
figure,
hold on,
grid on,
plot(lvec*L_application,S_avg(cnt),'-om');
title( ['Throughput VS Offered Load for ',num2str(N_stations),' Nodes']);
xlabel('Offered Load (bits/application/node)');
ylabel('Throughput (bits/application/node/s)');
figure,
hold on,
grid on,
plot(lvec*L_application,S_inst(cnt),'-og');
title( ['Instantaneous Throughput VS Offered Load for ',num2str(N_stations),' Nodes']);
xlabel('Offered Load (bits/application/node)');
ylabel('Instantaneous Throughput (bits/application/node/s)');
figure,
hold on,
grid on,
plot(lvec*L_application,Pcf(cnt),'-o');
hold on,
plot(lvec*L_application,Pcr(cnt),'-*');
title( ['Pcf/Pcr VS Offered Load for ',num2str(N_stations),' Nodes']);
legend('Pcf','Pcr');
xlabel('Offered Load (bits/application/node)');
ylabel('Pcf/Pcr');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of M-file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
