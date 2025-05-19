mutable struct MP_SV
    stable_bias::Float64
    mp_bias::Float64
    stable_cov::Float64
    mp_cov::Float64
    inovations::Circbuffer{Float64}
end


function mp_step!(mpsv::MP_SV, res, α)
    # mpsv - MP_SV object representing the multipath estimation filter of given satellite
    # res - residual of the measurement
    # α - threshold for the stable bias

    # returns: bias estimate, covariance of the estimate (and therefore the pseudorange measurement), status of the filter (1 = stable, 2 = updated bias, 3 = noisy)

    MPWINLEN = length(mpsv.inovations.data)

    push!(mpsv.inovations, res)
    stable_inovations = mpsv.inovations.data .- mpsv.stable_bias
    T_stable = sum(stable_inovations.^2) # sum of squares of inovations
    if T_stable < α * mpsv.stable_cov
        return mpsv.stable_bias, mpsv.stable_cov, 1
    end

    updated_stable_bias = ones(MPWINLEN)\mpsv.inovations.data
    updated_stable_inovations = mpsv.inovations.data .- updated_stable_bias
    T_update = sum(updated_stable_inovations.^2) # sum of squares of inovations
    if T_update < α * mpsv.stable_cov
        mpsv.stable_bias = updated_stable_bias
        return mpsv.stable_bias, mpsv.stable_cov, 2
    end
    σ2 = variance(mpsv.inovations.data)
    
    mpsv.mp_bias = ones(MPWINLEN) \ mpsv.inovations.data 
    mpsv.mp_cov = σ2

    meen = 1/2 * (mpsv.mp_bias + mpsv.stable_bias)
    kov = 1/2 * (mpsv.mp_cov + mpsv.stable_cov + (mpsv.mp_bias - meen)^2 + (mpsv.stable_bias - meen)^2)
    return meen, kov, 3
end


