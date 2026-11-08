# These tests moved from SMLMData 0.7.1's test/test_filters.jl. Since the move they compare positions in
# SMITE's pixels and expect the camera fields (XSize, YSize, PixelSize, NDims) in each saved file.
using SMLMio, Test
using SMLMCore: SMLMCore, IdealCamera, Emitter2DFit, Emitter3DFit, emitter_ndims
using MAT: matread, matwrite

mutable struct ForeignEmitter3DFit{T} <: SMLMCore.AbstractEmitter   # for save_smite
    x::T; y::T; z::T; photons::T; bg::T; σ_x::T; σ_y::T; σ_z::T; σ_photons::T; σ_bg::T
    frame::Int; dataset::Int; track_id::Int; id::Int
end
# Every field 0.7.0's save_smite reads, plus z, but no σ_z: 0.7.0 saved it without Z.
mutable struct ForeignZNoSigmaZ{T} <: SMLMCore.AbstractEmitter
    x::T; y::T; z::T; photons::T; bg::T; σ_x::T; σ_y::T; σ_photons::T; σ_bg::T
    frame::Int; dataset::Int; track_id::Int; id::Int
end
# Stores z but hides it: not in propertynames, refused by getproperty. Exposes σ_z. 0.7.0 saved it without Z.
mutable struct HiddenZLoc{T} <: SMLMCore.AbstractEmitter
    x::T; y::T; z::T; photons::T; bg::T; σ_x::T; σ_y::T; σ_z::T; σ_photons::T; σ_bg::T
    frame::Int; dataset::Int; track_id::Int; id::Int
end
Base.propertynames(::HiddenZLoc) = Tuple(f for f in fieldnames(HiddenZLoc) if f !== :z)
Base.getproperty(e::HiddenZLoc, s::Symbol) =
    s === :z ? throw(ArgumentError("HiddenZLoc keeps z private")) : getfield(e, s)
# Every field 0.7.0's save_smite reads, with z and σ_z of any type. MAT's writer throws on missing
# and Float16, so such a z must not reach the file; 0.7.0 saved these without Z.
struct OddZLoc{Z,S} <: SMLMCore.AbstractEmitter
    x::Float64; y::Float64; z::Z; photons::Float64; bg::Float64; σ_x::Float64; σ_y::Float64; σ_z::S
    σ_photons::Float64; σ_bg::Float64; frame::Int; dataset::Int; track_id::Int; id::Int
end
_oddz(z, σz) = OddZLoc(1.0, 1.5, z, 1000.0, 10.0, 0.01, 0.01, σz, 50.0, 2.0, 1, 1, 0, 1)

# Mimics SykTrack's Localization: position in a tuple, x/y/z computed in getproperty,
# :z listed in propertynames only for N >= 3. No emitter_ndims method of its own.
struct PropZLoc{N,T} <: SMLMCore.AbstractEmitter
    position::NTuple{N,T}
    uncertainty::NTuple{N,T}
    photons::T
    bg::T
    σ_photons::T
    σ_bg::T
    frame::Int
    dataset::Int
    track_id::Int
    id::Int
end
function Base.getproperty(e::PropZLoc{N,T}, s::Symbol) where {N,T}
    s === :x && return getfield(e, :position)[1]
    s === :y && return getfield(e, :position)[2]
    s === :z && return N >= 3 ? getfield(e, :position)[3] : zero(T)
    s === :σ_x && return getfield(e, :uncertainty)[1]
    s === :σ_y && return getfield(e, :uncertainty)[2]
    s === :σ_z && return N >= 3 ? getfield(e, :uncertainty)[3] : zero(T)
    return getfield(e, s)
end
Base.propertynames(::PropZLoc{N,T}) where {N,T} =
    (fieldnames(PropZLoc{N,T})..., (N >= 3 ? (:x, :y, :z, :σ_x, :σ_y, :σ_z) : (:x, :y, :σ_x, :σ_y))...)

# 0.7.0's save_smite body, copied, to compare the columns and values written.
function _save_smite_070(smld, filepath, filename)
    s = Dict{String,Any}()
    s["X"] = [e.x for e in smld.emitters]
    s["Y"] = [e.y for e in smld.emitters]
    if eltype(smld.emitters) <: Emitter3DFit
        s["Z"] = [e.z for e in smld.emitters]
    end
    s["Photons"] = [e.photons for e in smld.emitters]
    s["Bg"] = [e.bg for e in smld.emitters]
    s["X_SE"] = [e.σ_x for e in smld.emitters]
    s["Y_SE"] = [e.σ_y for e in smld.emitters]
    if eltype(smld.emitters) <: Emitter3DFit
        s["Z_SE"] = [e.σ_z for e in smld.emitters]
    end
    s["Photons_SE"] = [e.σ_photons for e in smld.emitters]
    s["Bg_SE"] = [e.σ_bg for e in smld.emitters]
    s["FrameNum"] = [e.frame for e in smld.emitters]
    s["DatasetNum"] = [e.dataset for e in smld.emitters]
    s["ConnectID"] = [e.track_id for e in smld.emitters]
    s["NFrames"] = smld.n_frames
    s["NDatasets"] = smld.n_datasets
    for (key, value) in smld.metadata
        haskey(s, key) || (s[key] = value)
    end
    matwrite(joinpath(filepath, filename), Dict("SMD" => s))
end

# Saves smld with 0.7.0's body and with save_smite (warnings silenced). Returns the keys the
# current save adds, or :changed if it drops or alters a 0.7.0 key, and the current file's SMD.
# A throw from save_smite propagates, failing the test. 0.7.0 wrote X, Y, X_SE and Y_SE in
# microns; save_smite writes SMITE's pixels, so those four compare after that conversion
# (cameras here start at 0).
function _legacy_compare(smld)
    mktempdir() do dir
        _save_smite_070(smld, dir, "old.mat")
        old = matread(joinpath(dir, "old.mat"))["SMD"]
        Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
            save_smite(smld, dir, "new.mat")
        end
        new = matread(joinpath(dir, "new.mat"))["SMD"]
        px = smld.camera.pixel_edges_x[2] - smld.camera.pixel_edges_x[1]
        pixels = Dict("X" => v -> v ./ px .+ 0.5, "Y" => v -> v ./ px .+ 0.5,
                      "X_SE" => v -> v ./ px, "Y_SE" => v -> v ./ px)
        same(k) = !haskey(pixels, k) ? isequal(new[k], old[k]) :
                  size(new[k]) == size(old[k]) && (isempty(old[k]) || new[k] ≈ pixels[k](old[k]))
        all(k -> haskey(new, k) && same(k), keys(old)) || return :changed, new
        return sort!(collect(setdiff(keys(new), keys(old)))), new
    end
end

@testset "Dimension routing" begin
    cam = IdealCamera(1:512, 1:512, 0.1)
    xs = [1.0, 5.0, 0.5, 1.0, 1.5]
    ys = [1.0, 1.0, 1.5, 9.0, 0.5]
    zs = [0.0, 0.0, 0.5, 0.0, 3.0]

    @testset "save_smite with foreign 3D emitters" begin
        es = [ForeignEmitter3DFit{Float64}(xs[i], ys[i], zs[i], 1000.0, 10.0, 0.01, 0.01, 0.02 + i,
                                           50.0, 2.0, i, 1, 0, i) for i in 1:5]
        s = SmiteSMLD{Float64,ForeignEmitter3DFit{Float64}}(es, cam, 5, 1, Dict{String,Any}())
        mktempdir() do dir
            save_smite(s, dir, "foreign.mat")
            smd = matread(joinpath(dir, "foreign.mat"))["SMD"]
            @test vec(smd["Z"]) == zs
            @test vec(smd["Z_SE"]) == [0.02 + i for i in 1:5]
        end
    end

    @testset "emitter_ndims: computed z via propertynames" begin
        cam = IdealCamera(1:64, 1:64, 0.1)
        mk2(x, y, i) = PropZLoc{2,Float64}((x, y), (0.01, 0.01), 100.0, 5.0, 1.0, 1.0, 1, 1, 0, i)
        mk3(x, y, z, i) = PropZLoc{3,Float64}((x, y, z), (0.01, 0.01, 0.02 + i), 100.0, 5.0, 1.0, 1.0, 1, 1, 0, i)
        pts = [(0.5, 0.5, 0.0), (1.5, 1.0, 0.5), (5.0, 1.0, 0.0), (1.0, 1.0, 3.0)]   # first two inside the box
        p2s = [mk2(p[1], p[2], i) for (i, p) in enumerate(pts)]
        p3s = [mk3(p..., i) for (i, p) in enumerate(pts)]

        mktempdir() do dir
            t3 = SmiteSMLD{Float64,PropZLoc{3,Float64}}(p3s, cam, 1, 1, Dict{String,Any}())
            save_smite(t3, dir, "p3.mat")
            smd = matread(joinpath(dir, "p3.mat"))["SMD"]
            @test vec(smd["Z"]) == [p[3] for p in pts]
            @test vec(smd["Z_SE"]) == [0.02 + i for i in 1:4]
            t2 = SmiteSMLD{Float64,PropZLoc{2,Float64}}(p2s, cam, 1, 1, Dict{String,Any}())
            save_smite(t2, dir, "p2.mat")
            smd = matread(joinpath(dir, "p2.mat"))["SMD"]
            @test !haskey(smd, "Z")
            @test !haskey(smd, "Z_SE")
        end
    end

    @testset "save_smite on mixed dimensions warns" begin
        cam = IdealCamera(1:64, 1:64, 0.1)
        A = SMLMCore.AbstractEmitter
        e2 = Emitter2DFit{Float64}(1.0, 1.0, 1000.0, 10.0, 0.01, 0.01, 50.0, 2.0)
        e3 = Emitter3DFit{Float64}(1.0, 1.0, 0.0, 1000.0, 10.0, 0.01, 0.01, 0.02, 50.0, 2.0)
        mixed = SmiteSMLD{Float64,A}(A[e2, e3], cam, 1, 1, Dict{String,Any}())
        empty_s = SmiteSMLD{Float64,A}(A[], cam, 1, 1, Dict{String,Any}())
        mktempdir() do dir
            # 0.7.0 saved mixed data without Z; it still does, and now says so.
            @test_logs (:warn, r"mixed 2D/3D") save_smite(mixed, dir, "mixed.mat")
            smd = matread(joinpath(dir, "mixed.mat"))["SMD"]
            @test vec(smd["X"]) ≈ [10.5, 10.5]   # x = 1.0 µm with 0.1 µm pixels: SMITE X = 10.5
            @test !haskey(smd, "Z")
            @test !haskey(smd, "Z_SE")
            @test_logs save_smite(empty_s, dir, "empty.mat")
            smd = matread(joinpath(dir, "empty.mat"))["SMD"]
            @test !haskey(smd, "Z")
        end
    end

    @testset "show labels mixed and unknown dimensions" begin
        cam = IdealCamera(1:64, 1:64, 0.1)
        A = SMLMCore.AbstractEmitter
        e2 = Emitter2DFit{Float64}(1.0, 1.0, 1000.0, 10.0, 0.01, 0.01, 50.0, 2.0)
        e3 = Emitter3DFit{Float64}(1.0, 1.0, 0.0, 1000.0, 10.0, 0.01, 0.01, 0.02, 50.0, 2.0)
        ms = SmiteSMLD{Float64,A}(A[e2, e3], cam, 1, 1, Dict{String,Any}())
        @test occursin("mixed 2D/3D", sprint(show, MIME("text/plain"), ms))
    end

    @testset "save_smite on an empty Union of fit types" begin
        # 0.7.0 wrote Z and Z_SE whenever the element type was a subtype of Emitter3DFit.
        cam = IdealCamera(1:64, 1:64, 0.1)
        U3 = Union{Emitter3DFit{Float32}, Emitter3DFit{Float64}}
        U2 = Union{Emitter2DFit{Float32}, Emitter2DFit{Float64}}
        mktempdir() do dir
            save_smite(SmiteSMLD{Float64,U3}(U3[], cam, 1, 1, Dict{String,Any}()), dir, "u3.mat")
            smd = matread(joinpath(dir, "u3.mat"))["SMD"]
            @test haskey(smd, "Z")
            @test haskey(smd, "Z_SE")
            save_smite(SmiteSMLD{Float64,U2}(U2[], cam, 1, 1, Dict{String,Any}()), dir, "u2.mat")
            smd = matread(joinpath(dir, "u2.mat"))["SMD"]
            @test !haskey(smd, "Z")
            @test !haskey(smd, "Z_SE")
        end
    end

    @testset "save_smite keeps 0.7.0's columns" begin
        cam = IdealCamera(1:64, 1:64, 0.1)
        A = SMLMCore.AbstractEmitter
        U3 = Union{Emitter3DFit{Float32}, Emitter3DFit{Float64}}
        e3 = Emitter3DFit{Float64}(1.0, 1.0, 0.5, 1000.0, 10.0, 0.01, 0.01, 0.125, 50.0, 2.0)
        e3b = Emitter3DFit{Float64}(2.0, 2.0, 0.75, 1000.0, 10.0, 0.01, 0.01, 0.5, 50.0, 2.0)
        e3f = Emitter3DFit{Float32}(1.0f0, 1.0f0, 0.25f0, 1000.0f0, 10.0f0, 0.01f0, 0.01f0, 0.0625f0,
                                    50.0f0, 2.0f0)
        fz = ForeignZNoSigmaZ{Float64}(1.0, 1.0, 0.5, 1000.0, 10.0, 0.01, 0.01, 50.0, 2.0, 1, 1, 0, 1)
        mktempdir() do dir
            # 0.7.0 dropped Z for an abstract-typed vector of 3D fits; now every emitter being
            # an Emitter3DFit writes it.
            sa = SmiteSMLD{Float64,A}(A[e3, e3b], cam, 1, 1, Dict{String,Any}())
            @test_logs save_smite(sa, dir, "a3.mat")
            smd = matread(joinpath(dir, "a3.mat"))["SMD"]
            @test vec(smd["Z"]) == [0.5, 0.75]
            @test vec(smd["Z_SE"]) == [0.125, 0.5]
            su = SmiteSMLD{Float64,U3}(U3[e3, e3f], cam, 1, 1, Dict{String,Any}())
            @test_logs save_smite(su, dir, "u3.mat")
            smd = matread(joinpath(dir, "u3.mat"))["SMD"]
            @test Float64.(vec(smd["Z"])) == [0.5, 0.25]
            @test Float64.(vec(smd["Z_SE"])) == [0.125, 0.0625]
            # A foreign emitter with z but no σ_z saves as in 0.7.0: no error, no Z, a warning.
            sf = SmiteSMLD{Float64,ForeignZNoSigmaZ{Float64}}([fz, fz], cam, 1, 1, Dict{String,Any}())
            @test_logs (:warn, r"3D emitters") save_smite(sf, dir, "fz.mat")
            smd = matread(joinpath(dir, "fz.mat"))["SMD"]
            @test vec(smd["X"]) ≈ [10.5, 10.5]   # x = 1.0 µm with 0.1 µm pixels: SMITE X = 10.5
            @test !haskey(smd, "Z")
            @test !haskey(smd, "Z_SE")
            # Mixed with 3D fits, or empty, it still saves as in 0.7.0: no Z.
            sm = SmiteSMLD{Float64,A}(A[e3, fz], cam, 1, 1, Dict{String,Any}())
            @test_logs (:warn, r"3D emitters") save_smite(sm, dir, "mz.mat")
            @test !haskey(matread(joinpath(dir, "mz.mat"))["SMD"], "Z")
            se = SmiteSMLD{Float64,ForeignZNoSigmaZ{Float64}}(ForeignZNoSigmaZ{Float64}[], cam, 1, 1,
                                                              Dict{String,Any}())
            @test_logs save_smite(se, dir, "ez.mat")
            @test !haskey(matread(joinpath(dir, "ez.mat"))["SMD"], "Z")
            # A stored z hidden from the properties, with σ_z exposed, is 2D and saves as in 0.7.0: no Z.
            hz = HiddenZLoc{Float64}(1.0, 1.0, 0.5, 1000.0, 10.0, 0.01, 0.01, 0.02, 50.0, 2.0, 1, 1, 0, 1)
            @test emitter_ndims(hz) == 2
            @test emitter_ndims([hz, hz]) == 2
            sh = SmiteSMLD{Float64,HiddenZLoc{Float64}}([hz, hz], cam, 1, 1, Dict{String,Any}())
            @test_logs save_smite(sh, dir, "hz.mat")
            smd = matread(joinpath(dir, "hz.mat"))["SMD"]
            @test vec(smd["X"]) ≈ [10.5, 10.5]   # x = 1.0 µm with 0.1 µm pixels: SMITE X = 10.5
            @test !haskey(smd, "Z")
            @test !haskey(smd, "Z_SE")
        end
    end

    @testset "save_smite against 0.7.0" begin
        # Where 0.7.0 saved, save_smite saves too, keeping every 0.7.0 key and value (positions in
        # SMITE's pixels); it adds the camera's NDims, PixelSize, XSize and YSize, and Z and Z_SE
        # for non-empty 3D data whose z and σ_z are all Float32 or Float64.
        camkeys = ["NDims", "PixelSize", "XSize", "YSize"]
        cam = IdealCamera(1:64, 1:64, 0.1)
        A = SMLMCore.AbstractEmitter
        sm(v) = SmiteSMLD{Float64,eltype(v)}(v, cam, 1, 1, Dict{String,Any}())
        e2 = Emitter2DFit{Float64}(1.0, 1.5, 1000.0, 10.0, 0.01, 0.01, 50.0, 2.0)
        e3 = Emitter3DFit{Float64}(1.0, 1.5, 0.5, 1000.0, 10.0, 0.01, 0.01, 0.125, 50.0, 2.0)
        e3f = Emitter3DFit{Float32}(1.0f0, 1.5f0, 0.25f0, 1000.0f0, 10.0f0, 0.01f0, 0.01f0, 0.0625f0,
                                    50.0f0, 2.0f0)
        f3 = ForeignEmitter3DFit{Float64}(1.0, 1.5, 0.5, 1000.0, 10.0, 0.01, 0.01, 0.02, 50.0, 2.0, 1, 1, 0, 1)
        f3b = ForeignEmitter3DFit{Float64}(2.0, 2.5, 0.75, 1000.0, 10.0, 0.01, 0.01, 0.04, 50.0, 2.0, 1, 1, 0, 2)
        e3b = Emitter3DFit{Float64}(2.0, 2.5, 0.75, 1000.0, 10.0, 0.01, 0.01, 0.25, 50.0, 2.0)
        fz = ForeignZNoSigmaZ{Float64}(1.0, 1.5, 0.5, 1000.0, 10.0, 0.01, 0.01, 50.0, 2.0, 1, 1, 0, 1)
        hz = HiddenZLoc{Float64}(1.0, 1.5, 0.5, 1000.0, 10.0, 0.01, 0.01, 0.02, 50.0, 2.0, 1, 1, 0, 1)
        U3 = Union{Emitter3DFit{Float32}, Emitter3DFit{Float64}}
        unchanged = (sm([e2, e2]), sm([e3, e3]), sm([e3f, e3f]), sm(U3[]), sm(A[]), sm(A[e2, e3]),
                     sm([fz, fz]), sm([hz, hz]), sm([_oddz(missing, 0.02), _oddz(missing, 0.02)]),
                     sm([_oddz(Float16(0.5), Float16(0.02)), _oddz(Float16(0.5), Float16(0.02))]),
                     sm(A[e3, _oddz(missing, 0.02)]),
                     # Float32 and Float64 z together make an AbstractFloat column, which MAT
                     # would write as a cell: no Z.
                     sm([_oddz(0.5, 0.02), _oddz(0.25f0, 0.03)]))
        for s in unchanged
            @test first(_legacy_compare(s)) == camkeys
        end
        widened = (sm(A[e3, e3b]), sm([f3, f3b]), sm([_oddz(0.5, 0.02), _oddz(0.25, 0.03)]),
                   sm([_oddz(0.5f0, 0.02f0), _oddz(0.25f0, 0.03f0)]))
        for s in widened
            added, new = _legacy_compare(s)
            @test added == [camkeys; "Z"; "Z_SE"]
            @test new["NDims"] == 3
            @test vec(new["Z"]) == [e.z for e in s.emitters]
            @test vec(new["Z_SE"]) == [e.σ_z for e in s.emitters]
            @test eltype(new["Z"]) === typeof(first(s.emitters).z)
        end
    end

    @testset "save_smite warns when data lose Z" begin
        cam = IdealCamera(1:64, 1:64, 0.1)
        A = SMLMCore.AbstractEmitter
        sm(v) = SmiteSMLD{Float64,eltype(v)}(v, cam, 1, 1, Dict{String,Any}())
        e2 = Emitter2DFit{Float64}(1.0, 1.5, 1000.0, 10.0, 0.01, 0.01, 50.0, 2.0)
        e3 = Emitter3DFit{Float64}(1.0, 1.5, 0.5, 1000.0, 10.0, 0.01, 0.01, 0.125, 50.0, 2.0)
        mktempdir() do dir
            @test_logs (:warn, r"mixed 2D/3D") save_smite(sm(A[e2, e3]), dir, "m.mat")
            @test_logs (:warn, r"3D emitters") save_smite(sm([_oddz(missing, 0.02)]), dir, "o.mat")
            @test_logs (:warn, r"3D emitters") save_smite(sm([_oddz(Float16(0.5), Float16(0.02))]), dir, "h.mat")
            @test_logs (:warn, r"3D emitters") save_smite(sm([_oddz(0.5, 0.02), _oddz(0.25f0, 0.03)]), dir, "f.mat")
            @test !haskey(matread(joinpath(dir, "h.mat"))["SMD"], "Z")
            @test_logs save_smite(sm([e2, e2]), dir, "2.mat")   # 2D data lose nothing
            @test_logs save_smite(sm(A[e3, e3]), dir, "3.mat")  # widened: Z written
            @test_logs save_smite(sm(A[]), dir, "e.mat")        # empty data lose nothing
        end
    end
end
