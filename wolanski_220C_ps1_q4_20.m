clear;
N = 500;
T = 20;
nreplic = 1000;
beta = 1;
beta_hat_20 = zeros(nreplic,1);
sigma_beta_hat_20 =  zeros(nreplic,1);
sigma_beta_twiddle_20 =  zeros(nreplic,1);

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
    beta_hat_20(j) = beta + left*right;
    uhat = y-x*beta_hat_20(j);
    uhatbar = zeros(1,N);
    for i = 1:N
        uhatbar(i) = mean(uhat(:,i));
    end
    demeaned_uhat = uhat - uhatbar;
    sigma_beta_twiddle_20(j) = (s_xx)^(-2)*(sum(sum(demeaned_x2.*((demeaned_uhat.^2)))));
    p = sum(demeaned_x.*demeaned_uhat);
    sigma_beta_hat_20(j) = (s_xx)^(-2)*sum(p.^2);
end

sd_beta_hat_20 = std(beta_hat_20);
sd_sigma_beta_hat_20 = std(sigma_beta_hat_20);
sd_sigma_beta_twiddle_20 = std(sigma_beta_twiddle_20);
E_sigma_beta_twiddle_20 = mean(sigma_beta_twiddle_20);
E_sigma_beta_hat_20 = mean(sigma_beta_hat_20);
bias_sigma_beta_hat_20 = E_sigma_beta_hat_20-sd_beta_hat_20;
bias_sigma_beta_twiddle_20 = E_sigma_beta_twiddle_20-sd_beta_hat_20;
rmse_sigma_beta_hat_20 = (mean(bias_sigma_beta_hat_20.^2))^(1/2);
rmse_sigma_beta_twiddle_20 = (mean(bias_sigma_beta_twiddle_20.^2))^(1/2);