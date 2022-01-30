using LsqFit

@doc """

	Parabolic WKB fit to the transformed electro-optic spectra
	F electric filed 
	L depletion width
	Eg bandgap

""" ->

q=1
m0=9.1e-31
m=0.276*m0
h=h/1.6e-19 # eV.s
ħ=h/(2*np.pi)


function WKB(En, yₕᵥₚ, p0=[1e+3, 0.01, 5.0] )

Yp(E,p) = sum(i->( abs.(q*p[3] .- q*E)./( ( ( (4*i^2 - 1)*( (ħ*(q*abs(p[1]) ).^i ) ).*(abs(p[2]*1e-7) ).^(i - 1) )/( sqrt(2*m)*2^(i + 1) ) ).^(2/(2*i + 1))  ) ).^(i + 1/2),1:2 )


#p0=[1e+3, 0.01, 5.0] # [F , L (100 nm), Eg(eV)]
fit=curve_fit(Yp,vec(En),vec(yₕᵥₚ), p0) 
cfp=coef(fit)
J=fit.jacobian
R=fit.resid

return vec(cf), J, R, vec(Yp(En,cfp))

end