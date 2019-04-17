function filter = MexiHatfilter(t,a,tau)

filter = 1/(sqrt(3*a)*pi^(1/4))*(1-((t-tau)/a).^2).*exp(1/(2*a^2)*-(t-tau).^2);

end

