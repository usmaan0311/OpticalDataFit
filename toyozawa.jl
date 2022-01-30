using LsqFit, QuadGK


@doc """
    Elliott-Toyozawa function used to fit absorption coeffiecient to extract excitonic
    binding energy and optical bandgap.

    Function is `toyozawa(En,α,d=10000,Rtol=1e-6,p0=[1,0.08,5.1,0.2,0.2])`

        En is energy vector
        α is absorption coefficient vector
        d is maximum discrete level to sum
        p0 is guess vector [c0 R Eg Γm Γc]
        c0 is normallization constant
        R is excitonic B.E
        Eg is bandgap
        Γm is Gaussian width of discrete excitonic spectra
        Γc is Gaussian width of excitonic continuum


        It returns a tuple of vec(cf), J, R, vec(αf(En,cf))

        cf is fitting coefficient vector
        J is Jacobian matrix
        R is Residual vector
        αf(En,cf) are fit values 
        
""" ->

function toyozawa(En,α,d=10000,Rtol=1e-6,p0=[1,0.08,5.1,0.2,0.2])

    αf(E,p) = ( p[1].*sqrt( abs(p[2]) )./(E.*(2π) )).*(sum(m-> (2*p[2]/(m^3))*(1/p[4])*exp.(-((E .- (p[3] - p[2]/(m^2))).^2)/(2*(p[4]^2))),1:d  ) .+ (1/p[5]).*quadgk(x-> exp.(-(x .- E).^2/(2*p[5].^2))./(1 .- exp.(-2π*sqrt( abs(p[2])./(abs(x .- p[3]) ))) ),p[3],Inf,rtol=Rtol)[1]  )
    
    αₙₒᵣₘ = α/maximum(α)
    
    #p0=[1,0.08,5.1,0.2,0.2] # [c0 R Eg Γm Γc]
    fit=curve_fit(αf,vec(En),αₙₒᵣₘ,p0)
    cf=coef(fit)
    J=fit.jacobian
    R=fit.resid
    
    return vec(cf), J, R, vec(αf(En,cf))*maximum(α)
        
end