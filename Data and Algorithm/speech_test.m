clear all
close all

load IR_sparse.mat;   % loads the echo path model
load speech.mat  % loads a speech signal

%% Genaration of synthetic speech sequence
 
% a = 0.8 ;
% Nr = [sqrt(1-a^2)] ;
% Dr = [1 -a] ;
% N = 50000 ;                 % sequence length
% u = randn(N,1);             
% far_end = filter(Nr,Dr,u);   % input sequence of unit variance (AR(1) process)
var_n = var(far_end)/100 ;               % Near end Background noise
N = length(far_end) ;

%% Impulse Response of the system

ho = IR_sparse ;
% ho = path ;

%% Echo/Desired signal

echo = filter(ho,1,far_end) ;
echo = echo + sqrt(var_n)*randn(N,1) ;

%% Initializations
M = length(ho) ;                    % Filter-Tap Length

disp('Please wait for a while...')
% IPNLMS Specifications
mu = 1 ;
alpha = 0 ;
delta = 9.8314e-04 ;
epsilon = 0.01 ;

w3 = zeros(M,1) ;            % Weight vector
u3 = zeros(M,1) ;            % regressor vector
m3 = zeros(N,1) ;            % misalignment vector


disp('Please wait for a while...')
%% IPNLMS Algorithm

for i=1:N
    u3 = [far_end(i);u3(1:M-1)];
    e3(i) = echo(i) - u3'*w3;
  
    for s = 1:M
        k(s) = (1-alpha)/2*M + (1+alpha)*norm(w3(s),1)/(2*norm(w3,1)+delta) ; % proportionate step-size implementation      
    end
    
    Q = diag(k) ;                  % Step-Size update matrix
    
    max_coeff2(i) = max(mu*k/(u3'*Q*u3 + epsilon)) ;
    min_coeff2(i) = min(mu*k/(u3'*Q*u3 + epsilon)) ;
    
    w3 = w3 + (mu*e3(i)*Q*u3)/(u3'*Q*u3 + epsilon) ;

    m3(i) = 20*log10(norm(ho-w3)/norm(ho)) ;
    
    if mod(i,5000)==0
         i
    end
end

% % OPLMS Specifications
% 
% w2 = zeros(M,1) ;            % Weight vector
% u2 = zeros(M,1) ;            % regressor vector
% m2 = zeros(N,1) ;            % misalignment matrix
% m = 1e-2 ;
% var_w = 0 ;                  % process noise
% gamma = ones(M,1);
% I_l = ones(M,1)  ;
% 
% %% OPLMS Algorithm
% 
% for i = 1:N
%    u2 = [far_end(i) ; u2(1:M-1)] ;
%    e2(i) = echo(i)-u2'*w2 ;
%    
%    var_x = (u2'*u2)/M  ;
%    q = M /(m + M*var_w);
%    mu_new = 1/(q*(var_n+var_x*(m+M*var_w))) ;
%    
%    w2_old = w2 ;
%    w2 = w2 + q*mu_new*(gamma.*u2)*e2(i) ;
%    
%    max_coeff1(i) = max(abs(q*mu_new*gamma)) ;
%    min_coeff1(i) = min(abs(q*mu_new*gamma)) ;
%    
%    m = m + M*var_w - q*mu_new*var_x*norm(gamma,2)^2 ;
%    
%    gamma = gamma + var_w*I_l + q*mu_new*var_x*(var_n+var_x*(m+M*var_w)-2*q*mu_new)*(gamma.*gamma) ;
%    r = max(q*gamma) ;
%    gamma = (1/r)*gamma ;
%    
%    var_w = (1/M)*norm(w2-w2_old,2)^2 ;
%   
%    
%    m2(i) = 20*log10(norm(ho-w2,2)/norm(ho,2)) ;
%     
%    if mod(i,5000)==0
%         i
%    end
% end



figure
subplot(311)
plot(1:N,far_end);
title('Far-end signal');
axis([0 N -2 2]);


subplot(312)
plot(1:N,echo);
title('Echo signal');
axis([0 N -2 2]);


subplot(313)
plot(1:N,e3);
title('Error signal');
axis([0 N -2 2]);
