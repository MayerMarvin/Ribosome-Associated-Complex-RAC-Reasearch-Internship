#!/usr/bin/env julia

module CoCoAssembly

using RCall # calculation of confidence intervals for HypothesisTests.FisherExactTest is _extremely_ slow

using Optim
using LineSearches
using Distributions
using Random
using Statistics

using DataStructures
using DataFrames

using Distributed
using HDF5
using ProgressMeter

export OptimDir,
       ll_max,
       ll_min,
       p_null_mle,
       p_null,
       p_ssig,
       p_dsig,
       ∇ssig,
       ∇dsig,
       ∇llssig,
       ∇llssig!,
       ∇lldsig,
       ∇lldsig!,
       binom_ll,
       TranscriptData,
       NullOptimizationResults,
       TranscriptResults,
       doFit,
       processData,
       postProcess

@enum OptimDir ll_max ll_min

function p_null_mle(mono::AbstractVector{<:Unsigned}, di::AbstractVector{<:Unsigned})
    sum(di) / (sum(mono) + sum(di))
end

function p_null(pars::AbstractVector{<:Real}, t::Real)
    @assert length(pars) == 1
    pars[1]
end

# pars[1]: Ibegin, pars[2]: Imax, pars[3]: a, pars[4]: tmid
function p_ssig_int(pars::AbstractVector{<:Real}, t::Real)
    @inbounds (pars[2] - pars[1]) / (1 + exp(-pars[3] * (t - pars[4]))) + pars[1]
end

function p_ssig(pars::AbstractVector{<:Real}, t::Real)
    @assert length(pars) == 4 "wrong dimension of parameter vector"
    p_ssig_int(pars, t)
end

# pars[1]: Ifinal, pars[2]: a, pars[3]: tmid, pars[4]: tdist
function p_ssig2_int(pars::AbstractVector{<:Real}, t::Real)
    @inbounds (1 - pars[1]) / (exp(-pars[2] * (t - (pars[3] + pars[4]))) + 1) + pars[1]
end

function p_ssig2(pars::AbstractVector{<:Real}, t::Real)
    @assert length(pars) == 4 "wrong dimension of parameter vector"
    p_ssig2_int(pars, t)
end

# pars[1]: Ibegin, pars[2]: Imax, pars[3]: Ifinal, pars[4]: a1, pars[5]: a2, pars[6]: tmid, pars[7]: tdist
function p_dsig(pars::AbstractVector{<:Real}, t::Real)
    @assert length(pars) == 7 "wrong dimension of parameter vector"
    @inbounds p_ssig_int(view(pars, [1, 2, 4, 6]), t) * p_ssig2_int(view(pars, [3, 5, 6, 7]), t)
end

# pars[1]: Ibegin, pars[2]: Imax, pars[3]: a, pars[4]: tmid
function ∇ssig(pars::AbstractVector{<:Real}, t::Real, dir::OptimDir=ll_max)
    @assert length(pars) == 4 "wrong dimension of parameter vector"
    @inbounds begin
        efun = exp(-pars[3] * (t - pars[4]))
        denom = 1 + efun

        ∂Imax = 1 / denom
        ∂Ibegin = 1 - 1 / denom
        ∂a = (pars[2] - pars[1]) * efun * (t - pars[4]) / denom^2
        ∂tmid = -(pars[2] - pars[1]) * efun * pars[3] / denom^2
        ∂t = (pars[2] - pars[1]) * efun * pars[3] / denom^2
        g = [∂Ibegin ∂Imax  ∂a ∂tmid ∂t]
        if (dir == ll_min)
            g .= -g
        end
        g
    end
end

function ∇ssig(pars::AbstractVector{<:Real}, t::AbstractVector{<:Real}, dir::OptimDir=ll_max)
    #hcat(∇ssig.((pars,), t)...)
    # benchmarked, this is twice as fast as the above version
    @assert length(pars) == 4 "wrong dimension of parameter vector"
    @inbounds begin
        efun = exp.(.-pars[3] .* (t .- pars[4]))
        denom = 1 .+ efun

        ∂Imax = 1 ./ denom
        ∂Ibegin = 1 .- ∂Imax
        ∂a = (pars[2] - pars[1]) .* efun .* (t .- pars[4]) ./ denom.^2
        ∂tmid = -(pars[2] - pars[1]) .* efun .* pars[3] ./ denom.^2
        ∂t = (pars[2] - pars[1]) .* efun .* pars[3] ./ denom.^2
        g = hcat(∂Ibegin, ∂Imax , ∂a, ∂tmid, ∂t)
        if dir == ll_min
            g .= -g
        end
        g
    end
end

# pars[1]: Ifinal, pars[2]: a, pars[3]: tmid, pars[4]: tdist
function ∇ssig2(pars::AbstractVector{<:Real}, t::Real)
    @assert length(pars) == 4 "wrong dimension of parameter vector"
    @inbounds begin
        efun = exp(-pars[2] * (t - (pars[3] + pars[4])))
        denom = efun + 1

        ∂Ifinal = 1 - 1 / denom
        ∂a = (1 - pars[1]) * efun * (t - (pars[3] + pars[4])) / denom^2
        ∂tmid = -(1 - pars[1]) * efun * pars[2] / denom^2
        ∂tdist = -(1 - pars[1]) * efun * pars[2] / denom^2
        ∂t = (1 - pars[1]) * efun * pars[2] / denom^2

        [∂Ifinal ∂a ∂tmid ∂tdist ∂t]
    end
end

function ∇ssig2(pars::AbstractVector{<:Real}, t::AbstractVector{<:Real})
    @assert length(pars) == 4 "wrong dimension of parameter vector"
    @inbounds begin
        efun = exp.(-pars[2] .* (t .- (pars[3] + pars[4])))
        denom = efun .+ 1

        ∂Ifinal = 1 .- 1 ./ denom
        ∂a = (1 - pars[1]) .* efun .* (t .- (pars[3] + pars[4])) ./ denom.^2
        ∂tmid = -(1 - pars[1]) .* efun .* pars[2] ./ denom.^2
        ∂tdist = -(1 - pars[1]) .* efun .* pars[2] ./ denom.^2
        ∂t = (1 - pars[1]) .* efun .* pars[2] ./ denom.^2

        hcat(∂Ifinal, ∂a, ∂tmid, ∂tdist, ∂t)
    end
end

# pars[1]: Ibegin, pars[2]: Imax, pars[3]: Ifinal, pars[4]: a1, pars[5]: a2, pars[6]: tmid, pars[7]: tdist
function ∇dsig(pars::AbstractVector{<:Real}, t, dir::OptimDir=ll_max)
    @assert length(pars) == 7 "wrong dimension of parameter vector"
    @inbounds begin
        ssig1 = p_ssig_int.((view(pars, [1, 2, 4, 6]),), t)
        ssig2 = p_ssig2_int.((view(pars, [3, 5, 6, 7]),), t)

        ∇1 = ∇ssig(view(pars, [1, 2, 4, 6]), t)
        ∇2 = ∇ssig2(view(pars, [3, 5, 6, 7]), t)

        ∂Ibegin = ∇1[:, 1] .* ssig2
        ∂Imax = ∇1[:, 2] .* ssig2
        ∂Ifinal = ∇2[:, 1] .* ssig1
        ∂a1 = ∇1[:, 3] .* ssig2
        ∂a2 = ∇2[:, 2] .* ssig1
        ∂tmid = ∇1[:, 4] .* ssig2 .+ ssig1 .* ∇2[:, 3]
        ∂tdist = ∇2[:, 4] .* ssig1
        ∂t = ∇1[:, 5] .* ssig2 .+ ssig1 .* ∇2[:, 5]

        g = hcat(∂Ibegin, ∂Imax, ∂Ifinal, ∂a1, ∂a2, ∂tmid, ∂tdist, ∂t)
        if dir == ll_min
            g .= -g
        end
        g
    end
end

# pars[1]: Ibegin, pars[2]: Imax, pars[3]: a, pars[4]: tmid
function ∇llssig(pars::AbstractVector{<:Real}, mono::AbstractVector{<:Unsigned}, di::AbstractVector{<:Unsigned}, t)
    @inbounds begin
        ssig = p_ssig.((pars,), t)
        ∇i = ∇ssig(pars, t)

        ∇o = di ./ ssig .- mono ./ (1 .- ssig)

        ∂Ibegin = ∇o .* ∇i[:, 1]
        ∂Imax = ∇o .* ∇i[:, 2]
        ∂a = ∇o .* ∇i[:, 3]
        ∂tmid = ∇o  .* ∇i[:, 4]
        [sum(∂Ibegin), sum(∂Imax), sum(∂a), sum(∂tmid)]
    end
end

function ∇llssig!(pars::AbstractVector{<:Real}, mono::AbstractVector{<:Unsigned}, di::AbstractVector{<:Unsigned}, t::AbstractVector{<:Real}, g::AbstractVector{<:Real}, dir::OptimDir=ll_max)
    @assert length(mono) == length(di) == length(t) "dimension mismatch"
    @assert length(pars) == 4 "wrong dimension of parameter vector"
    @assert length(g) == 4 "wrong dimension of gradient vector"
    g .= 0.
    @views @inbounds @simd for i in 1:length(t)
        efun = exp(-pars[3] * (t[i] - pars[4]))
        denom = 1 + efun
        ssig = (pars[2] - pars[1]) / denom + pars[1]#p_ssig_int(pars, t[i])

        ∇o = di[i] / ssig - mono[i] / (1 - ssig)

        g[1] += ∇o * (1 - 1 / denom)
        g[2] += ∇o / denom
        g[3] += ∇o * (pars[2] - pars[1]) * efun * (t[i] - pars[4]) / denom^2
        g[4] += ∇o * (-(pars[2] - pars[1]) * efun * pars[3] / denom^2)
    end
    if dir == ll_min
        g .= -g
    end
    nothing
end

# pars[1]: Ibegin, pars[2]: Imax, pars[3]: Ifinal, pars[4]: a1, pars[5]: a2, pars[6]: tmid, pars[7]: tdist
function ∇lldsig(pars::AbstractVector{<:Real}, mono::AbstractVector{<:Unsigned}, di::AbstractVector{<:Unsigned}, t, dir::OptimDir=ll_max)
    @inbounds begin
        p = p_dsig.((pars,), t)
        ∇i = ∇dsig(pars, t)

        ∇o = di ./ p .- mono ./ (1 .- p)

        ∂Ibegin = ∇o .* ∇i[:, 1]
        ∂Imax = ∇o .* ∇i[:, 2]
        ∂Ifinal = ∇o .* ∇i[:, 3]
        ∂a1 = ∇o .* ∇i[:, 4]
        ∂a2 = ∇o .* ∇i[:, 5]
        ∂tmid = ∇o .* ∇i[:, 6]
        ∂tdist = ∇o .* ∇i[:, 7]
        g = [sum(∂Ibegin), sum(∂Imax), sum(∂Ifinal), sum(∂a1), sum(∂a2), sum(∂tmid), sum(∂tdist)]
        if dir == ll_min
            g .= -g
        end
        g
    end
end

function ∇lldsig!(pars::AbstractVector{<:Real}, mono::AbstractVector{<:Unsigned}, di::AbstractVector{<:Unsigned}, t::AbstractVector{<:Real}, g::AbstractVector{<:Real}, dir::OptimDir=ll_max)
    @assert length(mono) == length(di) == length(t) "dimension mismatch"
    @assert length(pars) == 7 "wrong dimension of parameter vector"
    @assert length(g) == 7 "wrong dimension of gradient vector"
    g .= 0.
    @views @inbounds @simd for i in 1:length(t)
        efun1 = exp(-pars[4] * (t[i] - pars[6]))
        denom1 = 1 + efun1
        efun2 = exp(-pars[5] * (t[i] - (pars[6] + pars[7])))
        denom2 = efun2 + 1
        ssig1 = (pars[2] - pars[1]) / denom1 + pars[1]#p_ssig_int(pars[[1, 2, 4, 6]], t[i])
        ssig2 = (1 - pars[3]) / denom2 + pars[3]#p_ssig2_int(pars[[3, 5, 6, 7]], t[i])
        p = ssig1 * ssig2

        ∇o = di[i] / p - mono[i] / (1 - p[i])

        g[1] += ∇o * (1 - 1 / denom1) * ssig2
        g[2] += ∇o / denom1 * ssig2
        g[3] += ∇o * (1 - 1 / denom2) * ssig1
        g[4] += ∇o * (pars[2] - pars[1]) * efun1 * (t[i] - pars[6]) / denom1^2 * ssig2
        g[5] += ∇o * (1 - pars[3]) * efun2 * (t[i] - (pars[6] + pars[7])) / denom2^2 * ssig1
        g[6] += ∇o * ((-(pars[2] - pars[1]) * efun1 * pars[4] / denom1^2) * ssig2 + ssig1 * (-(1 - pars[3]) * efun2 * pars[5] / denom2^2))
        g[7] += ∇o * (-(1 - pars[3]) * efun2 * pars[5] / denom2^2) * ssig1
    end
    if dir == ll_min
        g .= -g
    end
    nothing
end

function binom_ll(pars::AbstractVector{<:Real}, mono::AbstractVector{<:Unsigned}, di::AbstractVector{<:Unsigned}, t::AbstractVector{<:Real}, pfun::Function, dir::OptimDir=ll_max)
    p = pfun.((pars,), t)
    if !all(zero.(p) .<= p .<= one.(p))
        if dir == ll_max
            return -Inf
        else
            return Inf
        end
    end
    ll = sum(logpdf.(Binomial.(di.+mono, p), di))
    if (dir == ll_min)
        ll = -ll
    end
    ll
end

struct TranscriptData
    pos::AbstractVector{<:Unsigned}
    mono::AbstractVector{<:Unsigned}
    di::AbstractVector{<:Unsigned}

    transcript::AbstractString
    gene_name::AbstractString
    cdsLength::Unsigned
    strand::AbstractString
end

mutable struct NullOptimizationResults{Tx, Tf} <: Optim.OptimizationResults
    minimizer::Tx
    minimum::Tf
end

mutable struct TranscriptResults
    nullfit::Optim.OptimizationResults
    ssigfit::Optim.OptimizationResults
    dsigfit::Optim.OptimizationResults

    bic::AbstractVector{<:Real}
    posterior_p::AbstractVector{<:Real}
    transcript::AbstractString
    gene::AbstractString
    cdsLength::Unsigned
    strand::AbstractString
    dataRange::Tuple{<:Unsigned, <:Unsigned}
    pos::AbstractVector{<:Unsigned}
    mono::AbstractVector{<:Unsigned}
    di::AbstractVector{<:Unsigned}
end

struct DatasetResults
    Σ_mono::UInt64
    Σ_di::UInt64
    fits::AbstractDict{<:AbstractString, TranscriptResults}
end

function randomInBound(bnd)
    (bnd[2] - bnd[1]) * rand() + bnd[1]
end

doFit(mono, di, t, pfun, gfun, bounds, defaultstart) = @views begin
    obj = pars -> binom_ll(pars, mono, di, t, pfun, ll_min)
    g! = (G, pars) -> gfun(pars, mono, di, t, G, ll_min)

    Random.seed!(42)
    fit = optimize(obj,
                g!,
                bounds[:, 1],
                bounds[:, 2],
                defaultstart,
                Fminbox(BFGS(linesearch=BackTracking())))
    while !Optim.converged(fit)
        fit = optimize(obj,
                g!,
                bounds[:, 1],
                bounds[:, 2],
                vec(mapslices(randomInBound, bounds, dims=2)),
                Fminbox(BFGS(linesearch=BackTracking())))
    end
    fit
end

##############################################################################
# postprocessing
##############################################################################

function Base.push!(d::DefaultDict, res::NullOptimizationResults)
    push!(d[:p], Optim.minimizer(res))
    push!(d[:loglik], -Optim.minimum(res))
end

function Base.push!(d::DefaultDict, res::Optim.MultivariateOptimizationResults, parnames::AbstractArray{<:AbstractString})
    pars = Optim.minimizer(res)
    @assert length(parnames) == length(pars)
    push!(d[:loglik], -Optim.minimum(res))
    for (n, p) in zip(parnames, pars)
        push!(d[Symbol(n)], p)
    end
end

function push_ftest(d::DefaultDict, res::TranscriptResults, ibefore::AbstractVector{<:Integer}, iafter::AbstractVector{<:Integer})
    dib = sum(res.di[ibefore])
    dia = sum(res.di[iafter])
    monob = sum(res.mono[ibefore])
    monoa = sum(res.mono[iafter])
    # TODO: replace with FisherExactTest from HypothesisTests.jl as soon as confint has reasonable speed
    ftest = rcopy(reval("fisher.test(matrix(c($dia, $dib, $monoa, $monob), ncol=2, byrow=TRUE), conf.level=0.95)"))
    push!(d[:or], ftest[:estimate])
    push!(d[:lo_CI], ftest[:conf_int][1])
    push!(d[:hi_CI], ftest[:conf_int][2])
end

function push_ssig(d::DefaultDict, res::TranscriptResults)
    push!(d, res.ssigfit, ["Ibegin", "Imax", "a", "tmid"])

    ibefore = findall(res.pos .< Optim.minimizer(res.ssigfit)[4])
    iafter = findall(res.pos .> Optim.minimizer(res.ssigfit)[4])
    push_ftest(d, res, ibefore, iafter)
end

function get_initial(guess::Number, bounds::Tuple{<:Number, <:Number})
    if guess - bounds[1] < eps()
        guess = bounds[1] + eps()
    elseif bounds[2] - guess < eps()
        guess = bounds[2] - eps()
    end
    guess
end

function push_dsig(d::DefaultDict, res::TranscriptResults)
    push!(d, res.dsigfit, ["Ibegin", "Imax", "Ifinal", "a1", "a2", "tmid", "tdist"])
    pars = Optim.minimizer(res.dsigfit)
    maxp = optimize(t -> -p_dsig(pars, t[1]), t -> [∇dsig(pars, t, ll_min)[1,8]], [Float64(res.dataRange[1])], [Float64(res.dataRange[2])], [Float64(res.dataRange[1] - 1 + findmax(p_dsig.((pars,), res.dataRange[1]:res.dataRange[2]))[2])], Fminbox(BFGS(linesearch=HagerZhang(linesearchmax=500))); inplace=false)
    push!(d[:maxt], Optim.minimizer(maxp)[1])
    push!(d[:maxp], -Optim.minimum(maxp))
    maxrng = (min(res.dataRange[1], 0.5 * pars[6]), Optim.minimizer(maxp)[1])
    minrng = (Optim.minimizer(maxp)[1], max(res.dataRange[2], 2 * (pars[6] + pars[7])))
    if Optim.minimizer(res.dsigfit)[1] > Optim.minimizer(res.dsigfit)[2]
        maxrng, minrng = minrng, maxrng
    end
    t1 = Optim.minimizer(optimize(t -> ∇dsig(pars, t, ll_min)[1, 8], [maxrng[1]], [maxrng[2]], [get_initial(pars[6], maxrng)], Fminbox(BFGS(linesearch=BackTracking())), autodiff=:forward))[1]
    t2 = Optim.minimizer(optimize(t -> ∇dsig(pars, t, ll_max)[1, 8], [minrng[1]], [minrng[2]], [get_initial(pars[6] + pars[7], minrng)], Fminbox(BFGS(linesearch=BackTracking())), autodiff=:forward))[1]
    if Optim.minimizer(res.dsigfit)[1] > Optim.minimizer(res.dsigfit)[2]
        t1, t2 = t2, t1
    end
    push!(d[:t1], t1)
    push!(d[:t2], t2)

    ibefore = findall(res.pos .< t1)
    iafter = findall(t1 .< res.pos .< t2)
    push_ftest(d, res, ibefore, iafter)
end

function push_fit(d::AbstractArray{<:DefaultDict}, fit::TranscriptResults, Σ_mono::Unsigned, Σ_di::Unsigned)
    @assert(length(d) == length(fit.bic))
    which = findmin(fit.bic)[2]
    if which == 1
        push!(d[1], fit.nullfit)
    elseif which == 2
        push_ssig(d[2], fit)
    else
        push_dsig(d[3], fit)
    end
    push!(d[which][:transcript], fit.transcript)
    push!(d[which][:gene], fit.gene)
    push!(d[which][:bic], fit.bic[which])
    push!(d[which][:posterior_p], fit.posterior_p[which])
    push!(d[which][:ndata], length(fit.pos))
    push!(d[which][:totalreads], sum(fit.mono) + sum(fit.di))
    push!(d[which][:totalRPM], sum(fit.mono) / Σ_mono * 1e6 + sum(fit.di) / Σ_di * 1e6)
    push!(d[which][:cdsLength], fit.cdsLength)
    push!(d[which][:minpos], fit.dataRange[1])
    push!(d[which][:maxpos], fit.dataRange[2])
    push!(d[which][:strand], fit.strand)
end

function DataFrames.DataFrame(d::DefaultDict{<:Any, <:Vector, <:Any})
    df = DataFrame()
    for k in sort(collect(keys(d)))
        df[k] = convert(Vector{typeof(d[k][1])}, d[k])
    end
    df
end

function normalize_ratio(r::Real, Σ_mono::Unsigned, Σ_di::Unsigned)
    r / (1 - r) * Σ_di / Σ_mono
end

function normalize_df!(df::DataFrame, Σ_mono::Unsigned, Σ_di::Unsigned, cols::AbstractVector{<:AbstractString})
    for par in cols
        df[Symbol(par * "_norm")] = normalize_ratio.(df[Symbol(par)], Σ_mono, Σ_di)
    end
    df
end

function normalize_ssig!(df::DataFrame, Σ_mono::Unsigned, Σ_di::Unsigned)
    normalize_df!(df, Σ_mono, Σ_di, ["Ibegin", "Imax", "a"])
end

function normalize_dsig!(df::DataFrame, Σ_mono::Unsigned, Σ_di::Unsigned)
    normalize_df!(df, Σ_mono, Σ_di, ["Ibegin", "Imax", "Ifinal", "a1", "a2"])
end

function postProcess(res::DatasetResults; progress::Bool=true)
    dicts = [DefaultDict{Symbol, Vector{Any}}(Vector{Any}) for i in 1:3]
    if progress
        pb = Progress(length(res.fits), desc="")
    end
    for (t, fit) in res.fits
        push_fit(dicts, fit, res.Σ_mono, res.Σ_di)
        progress && ProgressMeter.next!(pb; showvalues=[(:transcript, fit.transcript)])
    end
    base = DataFrame(dicts[1])
    ssig = DataFrame(dicts[2])
    dsig = DataFrame(dicts[3])
    (base=base, ssig=normalize_ssig!(ssig, res.Σ_mono, res.Σ_di), dsig=normalize_dsig!(dsig, res.Σ_mono, res.Σ_di))
end

##############################################################################
# multithreaded interface
##############################################################################

# Ibegin, Imax, a, tmid
function ssigbounds(t::Real)
    lb = [eps(), eps(), 0, 1]
    ub = [1 - eps(), 1 - eps(), 0.5, t]
    hcat(lb, ub)
end

# Ibegin, Imax, Ifinal, a1, a2, tmid, tdist
function dsigbounds(t::Real)
    lb = [eps(), eps(), eps(), 0, -0.5, 1, 1]
    ub = [1 - eps(), 1 - eps(), 1 - eps(), 0.5, 0, t, t]
    hcat(lb, ub)
end

function ssigbounds(t::AbstractVector{<:Real})
    ssigbounds(maximum(t))
end

function dsigbounds(t::AbstractVector{<:Real})
    dsigbounds(maximum(t))
end

function do_fits(in, out)
    try
        while (true)
            transcript = take!(in)
            ssigfit = doFit(transcript.mono,
                            transcript.di,
                            transcript.pos,
                            p_ssig,
                            ∇llssig!,
                            ssigbounds(transcript.pos),
                            [0.1, 0.9, 0.1, max(1, 0.5 * maximum(transcript.pos))])
            dsigfit = doFit(transcript.mono,
                            transcript.di,
                            transcript.pos,
                            p_dsig, ∇lldsig!,
                            dsigbounds(transcript.pos),
                            [0.1, 0.9, 0.1, 0.1, -0.1, max(1, 0.5 * maximum(transcript.pos)), max(1, 0.25 * maximum(transcript.pos))])

            pnull = p_null_mle(transcript.mono, transcript.di)
            nullmodel = binom_ll([pnull], transcript.mono, transcript.di, transcript.pos, p_null, ll_min)
            nullfit = NullOptimizationResults(pnull, nullmodel)

            ndata = length(transcript.pos)
            bics = [log(ndata) + 2 * Optim.minimum(nullfit),
                    log(ndata) * length(Optim.minimizer(ssigfit)) + 2 * Optim.minimum(ssigfit),
                    log(ndata) * length(Optim.minimizer(dsigfit)) + 2 * Optim.minimum(dsigfit)]
            posterior_p = exp.(-0.5 .* BigFloat.(bics))
            posterior_p = Float64.(posterior_p ./ sum(posterior_p))

            put!(out, (transcript.transcript, TranscriptResults(nullfit,
                        ssigfit,
                        dsigfit,
                        bics,
                        posterior_p,
                        transcript.transcript,
                        transcript.gene_name,
                        transcript.cdsLength,
                        transcript.strand,
                        extrema(transcript.pos),
                        transcript.pos,
                        transcript.mono,
                        transcript.di)))
        end
    catch y
        if isa(y, InterruptException)
            return
        else
            rethrow(y)
        end
    end
end

function setupWorkers(njobs=1)
    if length(workers()) > 1 || workers()[1] != 1
        rmprocs(workers())
    end
    if njobs > 1
        addprocs(njobs)
        @everywhere workers() include(@__FILE__)
    end

    jobs = RemoteChannel(() -> Channel{TranscriptData}(5*njobs))
    results = RemoteChannel(() -> Channel{Tuple{<:AbstractString, TranscriptResults}}(5*njobs))

    for p in workers()
        remote_do(do_fits, p, jobs, results)
    end
    (jobs=jobs, results=results)
end

function processData(infile::AbstractString, njobs::Unsigned=1; progress::Bool=true)
    jobs, results = setupWorkers(njobs)
    f = h5open(infile, "r")
    ndata = length(f)

    @async for t in names(f)
        dset=f[t]
        att = attrs(dset)
        data = read(dset)
        put!(jobs, TranscriptData(data[1,:],
                                  data[2,:],
                                  data[3,:],
                                  t,
                                  read(att, "gene_name"),
                                  read(att, "cds_length"),
                                  read(att, "strand")))
    end

    if progress
        pb = Progress(ndata, desc="")
    end

    res = Dict{String, TranscriptResults}()
    Σ_mono = 0
    Σ_di = 0
    for i in 1:ndata
        transcript, fit = take!(results)
        res[transcript] = fit
        Σ_mono += sum(fit.mono)
        Σ_di += sum(fit.di)
        progress && ProgressMeter.next!(pb; showvalues=[(:transcript, transcript)])
    end
    interrupt()
    DatasetResults(Σ_mono, Σ_di, res)
end

function processData(data::AbstractVector{TranscriptData}, njobs::Unsigned=1; progress::Bool=true)
    jobs, results = setupWorkers(njobs)
    @async for t in data
        put!(jobs, t)
    end
    ndata = length(data)
    if progress
        pb = Progress(ndata, desc="")
    end
    res = Dict{String, TranscriptResults}()
    Σ_mono = 0
    Σ_di = 0
    for i in 1:ndata
        transcript, fit = take!(results)
        res[transcript] = fit
        Σ_mono += sum(fit.mono)
        Σ_di += sum(fit.di)
        progress && ProgressMeter.next!(pb; showvalues=[(:transcript, transcript)])
    end
    interrupt()
    DatasetResults(Σ_mono, Σ_di, res)
end

end #module

using .CoCoAssembly
using Logging

using ArgParse
using FileIO
using Feather
if length(PROGRAM_FILE) > 0 && isfile(PROGRAM_FILE) && realpath(PROGRAM_FILE) == realpath(@__FILE__)
    global_logger(ConsoleLogger(stderr, Logging.Error))
    as = ArgParseSettings()
    @add_arg_table as begin
        "--procs", "-p"
            arg_type = Unsigned
            help = "number of processors to use"
            default = UInt(0)
        "INFILE"
            arg_type = AbstractString
            help = "path to input HDF5 file"
            required = true
    end
    args = parse_args(as, as_symbols=true)

    njobs = args[:procs]
    if njobs == 0
        njobs = UInt(Sys.CPU_THREADS * 0.5)
    end
    infile = args[:INFILE]

    inbase = splitext(infile)[1]
    outjld = inbase * ".jld2"

    res = processData(infile, njobs)

    save(outjld, basename(infile), res, compress=true)

    base, ssig, dsig = postProcess(res)
    Feather.write(inbase * "_base.feather", base)
    Feather.write(inbase * "_ssig.feather", ssig)
    Feather.write(inbase * "_dsig.feather", dsig)
end
