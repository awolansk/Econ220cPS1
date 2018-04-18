clear;
N = 500;
T = 5;
nreplic = 1000;
beta = 314.15926;
beta_hat = zeros(nreplic,1);
sigma_beta_hat =  zeros(nreplic,1);
sigma_beta_twiddle=  zeros(nreplic,1);

for j = 1:nreplic
    x = randn(T,N);
    u = zeros(T,N);
    for t = 1:T
        for i = 1:N
            u(t,i) = normrnd(0,(x(t,i))^2);
        end
    end
    y = x*beta+u;
    xbar = zeros(1,N);
    for i = 1:N
        xbar(i) = mean(x(:,i));
    end
    demeaned_x = x-xbar;
    demeaned_x2 = (demeaned_x).^2;
    s_xx = sum(sum(demeaned_x2));
    left = inv(s_xx);
    ubar = zeros(1,N);
    for i = 1:N
        ubar(i) = mean(u(:,i));
    end
    demeaned_u = u - ubar;
    demeaned_xu = demeaned_x .* demeaned_u;
    right = sum(sum(demeaned_xu));
    beta_hat(j) = beta + left*right;
    uhat = y-x*beta_hat(j);
    uhatbar = zeros(1,N);
    for i = 1:N
        uhatbar(i) = mean(uhat(:,i));
    end
    demeaned_uhat = uhat - uhatbar;
    sigma_beta_twiddle(j) = (s_xx)^(-2)*(sum(sum(demeaned_x2.*((demeaned_uhat.^2)))));
    p = sum(demeaned_x.*demeaned_uhat);
    sigma_beta_hat(j) = (s_xx)^(-2)*sum(p.^2);
end

sd_beta_hat = std(beta_hat);
sd_sigma_beta_hat = std(sigma_beta_hat);
sd_sigma_beta_twiddle = std(sigma_beta_twiddle);
E_sigma_beta_twiddle = mean(sigma_beta_twiddle);
E_sigma_beta_hat = mean(sigma_beta_hat);
bias_sigma_beta_hat = E_sigma_beta_hat-sd_beta_hat;
bias_sigma_beta_twiddle = E_sigma_beta_twiddle-sd_beta_hat;
rmse_sigma_beta_hat = (mean(bias_sigma_beta_hat.^2))^(1/2);
rmse_sigma_beta_twiddle = (mean(bias_sigma_beta_twiddle.^2))^(1/2);
