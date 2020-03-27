data { 
    int<lower=0> N;
    real t[N];
    real y[N];
    real<lower=0> sigma;
    real<lower=0, upper=pi()> theta1_prior_mean;
    real<lower=0, upper=pi()> theta2_prior_mean;
    real<lower=0> theta1_prior_sigma;
    real<lower=0> theta2_prior_sigma;
    real<lower=0, upper=50> f1_prior_mean;
    real<lower=300, upper=1000> f2_prior_mean;
    real<lower=0> f1_prior_sigma;
    real<lower=0> f2_prior_sigma;
    real<lower=0.5, upper=1> A1_prior_mean;
    real<lower=0, upper=0.1> A2_prior_mean;
    real<lower=0> A1_prior_sigma;
    real<lower=0> A2_prior_sigma;
}

parameters {
    real<lower=0, upper=pi()> theta1;
    real<lower=0, upper=pi()> theta2;
    real<lower=0> f1;
    real<lower=0> f2;
    real<lower=0> A1;
    real<lower=0> A2;
}

model {
    # generative model #
    for (n in 1:N)  {
        y[n]  ~  normal(A1*sin(2*pi()*f1*t[n] + theta1) + A2*sin(2*pi()*f2*t[n] + theta2), sigma);
    }
    # Prior #
    theta1 ~ normal(theta1_prior_mean, theta1_prior_sigma);
    theta2 ~ normal(theta2_prior_mean, theta2_prior_sigma);
    f1 ~ normal(f1_prior_mean, f1_prior_sigma);
    f2 ~ normal(f2_prior_mean, f2_prior_sigma);
    A1 ~ normal(A1_prior_mean, A1_prior_sigma);
    A2 ~ normal(A2_prior_mean, A2_prior_sigma);
}