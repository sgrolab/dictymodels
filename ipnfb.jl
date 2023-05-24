using StochasticDiffEq   # part of DifferentialEquations.jl

N = 100
Nt = 27
tspan = (0.0, 30*Nt)
cAMPin_level = 20
cAMPin_time = 30*Nt/2

ϵ = 0.1
γ = 0.5
c0 = 1.2
α = 0.058
α0 = 800.0
αPDE = 1000.0
Kd = 10^-5
S = 10^6
ρ = 10^-3.5
J = 0.5
σ = 0.15

u0 = [fill(-1.5,N); fill(-0.5,N); 0]

cAMPin(t) = t < cAMPin_time ? 0 : cAMPin_level  # step
#cAMPin(t) = t/tspan[2]*cAMPin_level  # ramp

I(x) = x/Kd>-1.0 ? α * log(1.0 + x/Kd) : NaN
Θ(x) = x>0.0 ? 1.0 : 0.0

iA = 1:N
iR = N+1:2N
icAMPe = 2N+1

function drift(du, u, p, t)
    @. du[iA] = u[iA] - u[iA]^3/3 - u[iR] + I(u[icAMPe])
    @. du[iR] = ϵ*(u[iA] - γ*u[iR] + c0)
    du[icAMPe] = cAMPin(t) + ρ*α0 + ρ*S/N*sum(Θ.(u[iA])) - (J + αPDE * ρ)*u[icAMPe]
end

function diffusion(du, u, p, t)
    du[iA] .= σ
    du[iR] .= 0
    du[icAMPe] = 0
end

prob = SDEProblem(drift, diffusion, u0, tspan)
sol = solve(prob, SOSRA())  # SOSRA is for problems with additive noise


using CairoMakie

fig = Figure()

ax = Axis(fig[1,1], ylabel = "A")
for i=1:10
    lines!(ax, sol.t./Nt, sol[iA[i],:])
end

ax = Axis(fig[2,1], ylabel = "R")
for i=1:10
    lines!(ax, sol.t./Nt, sol[iR[i],:])
end

ax = Axis(fig[3,1], ylabel = "cAMPe", xlabel = "time")
lines!(ax, sol.t./Nt, sol[icAMPe,:])

save("ipnfb-N100-rho$(ρ)-cAMPin$(cAMPin_level).pdf", fig)
