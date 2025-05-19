# DEPENDENCIES: 
# using LinearAlgebra

function timedate2unix(timedate::TimeDate)
    return datetime2unix(DateTime(timedate)) + Microsecond(timedate).value/1e6 + Nanosecond(timedate).value/1e9
end


function wls(A::AbstractArray{<:Real},y::Vector{<:Real},Σ::AbstractArray)
    # Weighted least squares
    # A: design matrix
    # y: observations
    # Σ: covariance matrix of observations
    # returns: x, covariance matrix of x
    # x = (A'Σ⁻¹A)⁻¹A'Σ⁻¹y
    # cov(x) = (A'Σ⁻¹A)⁻¹
    x = (A'*(Σ\A))\(A'/Σ)*y
    cov_x = (A'*(Σ\A))\I
    return x, cov_x
end


mutable struct Circbuffer{T}
    data::Vector{T}
    head::Int64
    size::Int64
end

function Circbuffer(size::Int64)
    return Circbuffer(zeros(size), 1, size)
end

function Base.push!(cb::Circbuffer, x::T) where T
    cb.data[cb.head] = x
    cb.head = mod1(cb.head + 1, cb.size)
end

function Base.getindex(cb::Circbuffer, i::Int64)
    return cb.data[mod1(cb.head - i, cb.size)]
end

function variance(x)
    return sum((x .- mean(x)).^2) / (length(x) - 1)
end

function ecef2enurotmat(lat, lon)
    lat = lat * π / 180.0
    lon = lon * π / 180.0
    return [
        -sin(lon) cos(lon) 0.0;
        -sin(lat)*cos(lon) -sin(lat)*sin(lon) cos(lat);
        cos(lat)*cos(lon) cos(lat)*sin(lon) sin(lat)
    ]
end

function enu2ecefrotmat(lat, lon)
    return ecef2enurotmat(lat, lon)'
end


function rms(x)
    return sqrt(mean(x.^2))
end
