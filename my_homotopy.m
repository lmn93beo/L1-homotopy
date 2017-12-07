% Some parameters
%rng(1234);
N = 512;   % signal length
M = round(N/2);    % no. of measurements
T = round(M/3);    % sparsity level
sType = 'randn'; % {'randn','sign','highD', 'blocks','pcwPoly'}
mType = 'randn';
SNR = 40;       % additive Gaussian noise

% Generate a signal
in = []; in.type = sType; in.T = T; in.randgen = 1;
x = genSignal(N,in); 

% measurement matrix
in = []; in.type = mType;
A = genAmat(M,N,in);

% measurements
sigma = sqrt(norm(A*x)^2/10^(SNR/10)/M);    

e = randn(M,1)*sigma;
y = A*x+e;

%% Homotopy
in.tau0 = max(abs(A' * y));
w = zeros(numel(x), 1);
J = find(abs(A' * y) == tau0);
in.maxiter = 10000;
in.tau_target = 0.03;
in.J = J;

out = do_homotopy(A, y, in);


% for i = 1:maxiter
%     JC = 1 : numel(x);
%     JC(J) = [];
% 
%     % Compute delta x
%     eta = sign(A' * (y - A * w));
%     AJ = A(:,J);
%     AJC = A(:,JC);
%     dw = zeros(numel(x), 1);
%     dw(J) = (AJ' * AJ)^(-1) * eta(J);
% 
%     % Smallest step to change in support
%     stepJ = -w(J) ./ dw(J);
%     p = AJC' * (A * w - y);
%     q = AJC' * (A * dw);
%     step_JC_plus = (tau0 - p) ./ (1 + q);
%     step_JC_minus = (-tau0 - p) ./ (q - 1);
% 
%     step_all = [stepJ; step_JC_plus; step_JC_minus];
%     step_all = step_all(step_all > 0);
%     min_step = min(step_all);
% 
%     idplus = find(step_JC_plus == min_step);
%     idminus = find(step_JC_minus == min_step);
%     idJ = find(stepJ == min_step);
% 
%     if ~isempty(idplus)
%         %fprintf('Added entry %d\n', JC(idplus));
%         J = [J JC(idplus)];
%     elseif ~isempty(idminus)
%         %fprintf('Added entry %d\n', JC(idminus));
%         J = [J JC(idminus)];
%     else
%         %fprintf('Removed entry %d\n', J(idJ));
%         J(idJ) = [];
%     end
%     %disp(J);
% 
%     tau0 = tau0 - min_step;
%     w = w + min_step * dw;
%     
%     if tau0 < tau_target
%         break;
%     end
%     
% end
% 
plot(x);
hold on;
plot(out.w);