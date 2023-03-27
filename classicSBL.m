function [errorc,time,mu_xc_est,gammare,H_re] = classicSBL(numItr,H,Res1,noise_var,y,K2,irs_pattern,A_irs_a,A_irs_d,A1,A2,H1,H2,SNR,thres)

time = 0;

gammac = 1*ones(Res1^3,1); % initialize the prior
mu_x_c = 1*ones(Res1^3,1);
keep_listc = 1:Res1^3;

normc = 1;
itrc = 1;


if SNR > 13 && SNR < 25
    r = 1e-4;
elseif SNR >= 25
    r = 1e-3;
else
    r = 1e-4;
end

while(itrc < numItr + 1 && normc > thres) % do the iteration
    itrc
    
    tic;
    
    gammac_t = zeros(Res1^3,1);
    gammac_t(keep_listc) = gammac;
    gammac_old = gammac_t;

    mu_xc_t = zeros(Res1^3,1);
    mu_xc_t(keep_listc) = mu_x_c;
    mu_xc_old = mu_xc_t;

    if min(abs(gammac)) < r*max(abs(gammac)) % 
         
        indexc = find(gammac > r*max(abs(gammac))); % should keep
        gammac = gammac(indexc);
        H = H(:,indexc);
        keep_listc = keep_listc(indexc);   
    end

    Gammac = diag(gammac);
    G = diag(sqrt(gammac));
    % update the posterior mean and variance
    
    invPhic = H'* (noise_var * eye(size(H,1)) + H * Gammac * H')^-1 * H;
    Sigma_xc = Gammac - Gammac * invPhic * Gammac;
    mu_x_c = noise_var^-1 * Sigma_xc * H' * y;

    
    % updating
    Sigma_xc_diag = real(diag(Sigma_xc));

    gammac = Sigma_xc_diag + abs(mu_x_c).^2;
    Gammac = diag(gammac);
    G = diag(sqrt(gammac));
   
    time = time + toc;
    % re-estimating for error computation
    
    invPhic = H'* (noise_var * eye(size(H,1)) + H * Gammac * H')^-1 * H;
    Sigma_xc = Gammac - Gammac * invPhic * Gammac;
    mu_x_c = noise_var^-1 * Sigma_xc * H' * y;

    mu_xc_est = zeros(Res1^3,1);
    mu_xc_est(keep_listc) = mu_x_c;
    gammare = zeros(Res1^3,1);
    gammare(keep_listc) = gammac;

    normc = norm(mu_xc_est-mu_xc_old)/norm(mu_xc_old)

    % computing the error for each IRS pattern and averaging them all
    for i = 1:K2
        Htt = irs_pattern(:,i).'*kr(A_irs_a.',A_irs_d').';
        Ht(:,:,i) = Htt(:,1:Res1);
        H_re(:,:,i) = kron(kron(Ht(:,:,i),conj(A1)),A2)*mu_xc_est;
        Htrue(:,:,i) = vec(H2*diag(irs_pattern(:,i))*H1);       
        nmse(itrc,i) = norm(H_re(:,:,i) - Htrue(:,:,i),'fro')/norm(Htrue(:,:,i),'fro');
    end
    
    error = sum(nmse(itrc,:))/K2
    errorc(itrc) = error;

    itrc = itrc + 1; 
end
errorc = errorc(end);
end

