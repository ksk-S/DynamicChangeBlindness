data { 
    int<lower=0> N;
    real t[N];
    real y[N];
    real<lower=0> sigma;
    real<lower=-1.0, upper=1.0> theta1_prior_mean;
    real<lower=-1.0, upper=1.0> theta2_prior_mean;
    real<lower=0> theta1_prior_sigma;
    real<lower=0> theta2_prior_sigma;
    real<lower=0.0, upper=1.0> f2_prior_mean;
    real<lower=0> f2_prior_sigma;
    real<lower=0, upper=30> f1;
    real<lower=0> A1;
    real<lower=0> A2;
}

parameters {
    real<lower=-1.0, upper=1.0> theta1;
    real<lower=-1.0, upper=1.0> theta2;
    real<lower=0.0, upper=1.0> f2;
}

transformed parameters {
    real dt;
     dt = t[2] - t[1];
}

model {
    # generative model #
    for (n in 1:N)  {
        y[n]  ~  normal(A1*sin(2*pi()*f1*t[n] + theta1*pi()) + A2*sin(2*pi()*(f2*200)*t[n] + theta2*pi()), sigma);
    }
    # Prior #
    theta1 ~ normal(theta1_prior_mean, theta1_prior_sigma);
    theta2 ~ normal(theta2_prior_mean, theta2_prior_sigma);
    f2 ~ normal(f2_prior_mean, f2_prior_sigma);
}

generated quantities {
    real y_pred[N];
    for (n in 1:N) {
        y_pred[n] = A1*sin(2*pi()*f1*t[n] + theta1*pi()) + A2*sin(2*pi()*(f2*500)*t[n] + theta2*pi());
    }
}