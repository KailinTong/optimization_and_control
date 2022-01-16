function value = z_pi(alpha_i, beta_i, gamma_i, u_i)
if (beta_i + gamma_i*u_i) > 0
    value = alpha_i + sqrt(beta_i + gamma_i*u_i);
else
    value = 0;
end

end