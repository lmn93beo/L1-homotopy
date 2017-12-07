function out = do_homotopy(A, y, in)
tau0 = in.tau0;
tau_target = in.tau_target;
maxiter = in.maxiter;
xlen = size(A, 2);
w = zeros(xlen, 1);
J = find(abs(A' * y) == tau0);


for i = 1:maxiter
    JC = 1 : xlen;
    JC(J) = [];

    % Compute delta x
    eta = sign(A' * (y - A * w));
    AJ = A(:,J);
    AJC = A(:,JC);
    dw = zeros(xlen, 1);
    dw(J) = (AJ' * AJ)^(-1) * eta(J);

    % Smallest step to change in support
    stepJ = -w(J) ./ dw(J);
    p = AJC' * (A * w - y);
    q = AJC' * (A * dw);
    step_JC_plus = (tau0 - p) ./ (1 + q);
    step_JC_minus = (-tau0 - p) ./ (q - 1);

    step_all = [stepJ; step_JC_plus; step_JC_minus];
    step_all = step_all(step_all > 0);
    min_step = min(step_all);

    idplus = find(step_JC_plus == min_step);
    idminus = find(step_JC_minus == min_step);
    idJ = find(stepJ == min_step);
    
%     J = [J JC(step_JC_plus == min_step)];
%     J = [J JC(step_JC_minus == min_step)];
%     J = [J JC(stepJ == min_step)];

    if ~isempty(idplus)
        %fprintf('Added entry %d\n', JC(idplus));
        J = [J JC(idplus)];
    elseif ~isempty(idminus)
        %fprintf('Added entry %d\n', JC(idminus));
        J = [J JC(idminus)];
    else
        %fprintf('Removed entry %d\n', J(idJ));
        J(idJ) = [];
    end
    %disp(J);

    tau0 = tau0 - min_step;
    w = w + min_step * dw;
    
    if tau0 < tau_target
        break;
    end
    
end

out.w = w;
out.tau0 = tau0;
out.iters = i;