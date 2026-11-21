# These tests moved from SMLMData 0.7.1's test/test_filters.jl. Since the move they compare positions in
# SMITE's pixels and expect the camera fields (XSize, YSize, PixelSize, NDims) in each saved file.
using SMLMio, Test
using SMLMCore: SMLMCore, IdealCamera, Emitter2DFit, Emitter3DFit, emitter_ndims
using MAT: matread, matwrite, matopen

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

# A 2D fit emitter whose fields may be any number type, and an SMLD whose counts may be.
struct AnyNumLoc <: SMLMCore.AbstractEmitter
    x::Any; y::Any; photons::Any; bg::Any; σ_x::Any; σ_y::Any; σ_photons::Any; σ_bg::Any
    frame::Any; dataset::Any; track_id::Any; id::Int
end
_anyloc(; x=1.0, photons=1000.0, bg=10.0, frame=1, dataset=1, track_id=0) =
    AnyNumLoc(x, 1.5, photons, bg, 0.01, 0.01, 50.0, 2.0, frame, dataset, track_id, 1)
struct AnyCountSMLD{E} <: SMLMCore.AbstractSMLD
    emitters::Vector{E}; camera::SMLMCore.AbstractCamera; n_frames::Any; n_datasets::Any
    metadata::Dict{String,Any}
end

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
            # Float32 and Float64 emitters together: every column is a double array, never the
            # cell MAT writes for a mixed vector, and the file loads back through the loader.
            cols = ["X", "Y", "Z", "X_SE", "Y_SE", "Z_SE", "Photons", "Bg", "Photons_SE", "Bg_SE",
                    "FrameNum", "DatasetNum", "ConnectID"]
            @test all(k -> smd[k] isa AbstractArray{Float64}, cols)
            back = load_smite_3d(SmiteSMD(dir, "u3.mat"))
            @test [e.z for e in back.emitters] == [0.5, 0.25]
            @test [e.σ_z for e in back.emitters] == [0.125, 0.0625]
            @test [e.photons for e in back.emitters] == [1000.0, 1000.0]
            @test [e.x for e in back.emitters] ≈ [1.0, 1.0]
            # The same in 2D, and with a Float32 camera, whose pixel makes the position columns mixed too
            U2 = Union{Emitter2DFit{Float32}, Emitter2DFit{Float64}}
            e2 = Emitter2DFit{Float64}(1.0, 2.0, 900.0, 9.0, 0.01, 0.02, 40.0, 1.0)
            e2f = Emitter2DFit{Float32}(3.0f0, 4.0f0, 800.0f0, 8.0f0, 0.01f0, 0.02f0, 30.0f0, 1.5f0)
            for c in (cam, IdealCamera(1:64, 1:64, 0.1f0))
                save_smite(SmiteSMLD{Float64,U2}(U2[e2, e2f], c, 1, 1, Dict{String,Any}()), dir, "u2.mat")
                smd = matread(joinpath(dir, "u2.mat"))["SMD"]
                @test all(k -> smd[k] isa AbstractArray{Float64}, setdiff(cols, ["Z", "Z_SE"]))
                back = load_smite_2d(SmiteSMD(dir, "u2.mat"))
                @test [e.photons for e in back.emitters] == [900.0, 800.0]
                @test [e.σ_bg for e in back.emitters] == [1.0, 1.5]
                @test [e.y for e in back.emitters] ≈ [2.0, 4.0] rtol = 1e-6
            end
            # Only Float32 values in an abstract-typed vector stay single, as Float32 data do
            save_smite(SmiteSMLD{Float64,U2}(U2[e2f, e2f], cam, 1, 1, Dict{String,Any}()), dir, "u2f.mat")
            smd = matread(joinpath(dir, "u2f.mat"))["SMD"]
            @test all(k -> eltype(smd[k]) == Float32, ["X", "Photons", "Bg", "X_SE", "Photons_SE", "Bg_SE"])
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
        # for non-empty 3D data whose z and σ_z are real numbers a double holds exactly (single
        # when all Float32, otherwise double).
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
                     sm(A[e3, _oddz(missing, 0.02)]))
        for s in unchanged
            @test first(_legacy_compare(s)) == camkeys
        end
        # Float32 and Float64 z together make an AbstractFloat column, which MAT would write as a
        # cell, and MAT cannot write Float16 at all; both are written as double.
        widened = (sm(A[e3, e3b]), sm([f3, f3b]), sm([_oddz(0.5, 0.02), _oddz(0.25, 0.03)]),
                   sm([_oddz(0.5f0, 0.02f0), _oddz(0.25f0, 0.03f0)]),
                   sm([_oddz(0.5, 0.02), _oddz(0.25f0, 0.03)]),
                   sm([_oddz(Float16(0.5), Float16(0.02)), _oddz(Float16(0.25), Float16(0.02))]))
        for s in widened
            added, new = _legacy_compare(s)
            @test added == [camkeys; "Z"; "Z_SE"]
            @test new["NDims"] == 3
            @test vec(new["Z"]) == [e.z for e in s.emitters]
            @test vec(new["Z_SE"]) == [e.σ_z for e in s.emitters]
            @test eltype(new["Z"]) === (all(e -> e.z isa Float32, s.emitters) ? Float32 : Float64)
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
            # z that are not real numbers save without Z (missing, above); a real z a double would
            # round (BigFloat) is an error naming the emitter, and nothing is written
            @test_throws r"z of emitter 2 is .*MATLAB double holds exactly" save_smite(
                sm([_oddz(0.5, 0.02), _oddz(big"0.1", 0.02)]), dir, "b.mat")
            @test !isfile(joinpath(dir, "b.mat"))
            @test_throws r"σ_z of emitter 1" save_smite(sm([_oddz(0.5, big"0.1")]), dir, "b.mat")
            # A complex z with no imaginary part is its real part, as 0.1.0 saved it
            @test_logs save_smite(sm([_oddz(0.25 + 0.0im, 0.02)]), dir, "c.mat")
            @test matread(joinpath(dir, "c.mat"))["SMD"]["Z"] === 0.25
            # Float16 z loses nothing: written as double
            @test_logs save_smite(sm([_oddz(Float16(0.5), Float16(0.02))]), dir, "h.mat")
            @test matread(joinpath(dir, "h.mat"))["SMD"]["Z"] === 0.5
            # Float32 and Float64 z together lose nothing: written as double
            @test_logs save_smite(sm([_oddz(0.5, 0.02), _oddz(0.25f0, 0.03)]), dir, "f.mat")
            @test eltype(matread(joinpath(dir, "f.mat"))["SMD"]["Z"]) == Float64
            @test [e.z for e in load_smite_3d(SmiteSMD(dir, "f.mat")).emitters] == [0.5, 0.25]
            @test_logs save_smite(sm([e2, e2]), dir, "2.mat")   # 2D data lose nothing
            @test_logs save_smite(sm(A[e3, e3]), dir, "3.mat")  # widened: Z written
            @test_logs save_smite(sm(A[]), dir, "e.mat")        # empty data lose nothing
        end
    end
end

@testset "save_smite checks keys first and never destroys a file" begin
    cam = IdealCamera(1:64, 1:64, 0.1)
    e = Emitter2DFit{Float64}(1.0, 1.5, 1000.0, 10.0, 0.01, 0.01, 50.0, 2.0)
    mk(md) = SmiteSMLD{Float64,Emitter2DFit{Float64}}([e], cam, 1, 1, md)
    good = mk(Dict{String,Any}("analysis_date" => "today"))
    mktempdir() do dir
        save_smite(good, dir, "r.mat")
        before = read(joinpath(dir, "r.mat"))

        @testset "invalid keys are named, one error, file untouched" begin
            md = Dict{String,Any}("analysis date" => 1.0, "γ" => 2.0, "1st_pass" => 3.0,
                                  "a" ^ 64 => 4.0, "drift-model" => 5.0, "fine_key" => 6.0,
                                  "monomer_D" => Dict(1 => 0.1, 2 => 0.2),
                                  "nested" => Dict("ok" => 1.0, "bad key" => 2.0))
            err = try save_smite(mk(md), dir, "r.mat"); nothing catch ex; ex end
            @test err isa ArgumentError
            msg = err.msg
            for k in ("analysis date", "γ", "1st_pass", "a" ^ 64, "drift-model", "monomer_D", "bad key")
                @test occursin(k, msg)
            end
            @test occursin("key 1 of the Dict monomer_D", msg)
            @test occursin("key 2 of the Dict monomer_D", msg)
            @test !occursin("fine_key", msg)
            @test read(joinpath(dir, "r.mat")) == before
            @test readdir(dir) == ["r.mat"]
            # No file is created either when there was none
            @test_throws ArgumentError save_smite(mk(md), dir, "new.mat")
            @test !isfile(joinpath(dir, "new.mat"))
            @test readdir(dir) == ["r.mat"]
        end

        @testset "a Dict with Symbol keys is named too" begin
            @test_throws r"sym" save_smite(mk(Dict{String,Any}("sym" => Dict(:a => 1.0))), dir, "r.mat")
            @test read(joinpath(dir, "r.mat")) == before
        end

        @testset "a write that fails midway leaves the existing file" begin
            # MAT cannot write Float16: the failure comes from matwrite, after the keys check
            @test_throws Exception save_smite(mk(Dict{String,Any}("h" => Float16(1))), dir, "r.mat")
            @test read(joinpath(dir, "r.mat")) == before
            @test readdir(dir) == ["r.mat"]
        end

        @testset "normal saves are unaffected" begin
            ok = mk(Dict{String,Any}("analysis_date" => "later", "a" => "a" ^ 100,
                                     "z" ^ 63 => 1.0, "inner" => Dict("k_1" => 2.0)))
            @test_logs save_smite(ok, dir, "r.mat")   # overwrites the existing file
            smd = matread(joinpath(dir, "r.mat"))["SMD"]
            @test smd["analysis_date"] == "later"
            @test smd["inner"]["k_1"] == 2.0
            @test haskey(smd, "z" ^ 63)
            @test load_smite_2d(SmiteSMD(dir, "r.mat")).emitters[1].x ≈ 1.0
            @test readdir(dir) == ["r.mat"]
        end

        @testset "a key with a trailing newline; arrays of Dicts" begin
            now = read(joinpath(dir, "r.mat"))
            @test_throws r"a\\n" save_smite(mk(Dict{String,Any}("a\n" => 1.0)), dir, "r.mat")
            @test read(joinpath(dir, "r.mat")) == now
            # MAT writes an array of Dicts as a struct array: its keys must be MATLAB names and
            # strings (MAT 0.10 refuses Symbol keys there; newer MAT takes them)
            @test_throws r"bad key" save_smite(mk(Dict{String,Any}("arr" => [Dict("bad key" => 1.0)])), dir, "r.mat")
            @test_throws r":a" save_smite(mk(Dict{String,Any}("arr" => [Dict(:a => 1.0)])), dir, "r.mat")
            @test read(joinpath(dir, "r.mat")) == now
            save_smite(mk(Dict{String,Any}("arr" => [Dict("a" => 1.0), Dict("a" => 2.0)])), dir, "r.mat")
            @test haskey(matread(joinpath(dir, "r.mat"))["SMD"], "arr")
            @test readdir(dir) == ["r.mat"]
        end

        # The temporary copies live in the system's temporary folder: none may be left there
        ours() = filter(startswith("smlmio-$(getpid())-"), readdir(tempdir()))
        @testset "a new file whose write fails is not left behind" begin
            @test_throws Exception save_smite(mk(Dict{String,Any}("h" => Float16(1))), dir, "new.mat")
            @test readdir(dir) == ["r.mat"]
            @test isempty(ours())
        end

        @testset "a file held open by HDF5 is refused and left as it was" begin
            # A read-only handle: closing a read-write one updates the file's HDF5 header itself
            now = read(joinpath(dir, "r.mat"))
            io = matopen(joinpath(dir, "r.mat"), "r")
            try
                @test_throws Exception save_smite(good, dir, "r.mat")
            finally
                close(io)
            end
            @test read(joinpath(dir, "r.mat")) == now
            @test load_smite_2d(SmiteSMD(dir, "r.mat")).emitters[1].x ≈ 1.0
            @test isempty(ours())
        end

        @testset "too little free space is refused before the file is touched" begin
            now = read(joinpath(dir, "r.mat"))
            @test_throws r"more bytes" SMLMio._check_space(joinpath(dir, "r.mat"), typemax(Int))
            @test read(joinpath(dir, "r.mat")) == now
        end

        @testset "a final copy that fails: a new file is removed, no temporary file stays" begin
            failing(src, dst) = (write(dst, "partial"); error("disk full"))
            s = Dict{String,Any}("X" => [1.0])
            @test_throws r"disk full" SMLMio._write_smd(joinpath(dir, "new2.mat"), s; copy=failing)
            @test !ispath(joinpath(dir, "new2.mat"))
            @test isempty(ours())
            # An existing file is then partly written, as documented, and the error says so; the
            # temporary file goes
            err = try SMLMio._write_smd(joinpath(dir, "r.mat"), s; copy=failing); nothing catch ex; ex end
            @test err isa ErrorException
            @test occursin("incomplete", err.msg) && occursin("save again", err.msg) &&
                  occursin("disk full", err.msg)
            @test read(joinpath(dir, "r.mat"), String) == "partial"
            @test isempty(ours())
        end
    end
end

@testset "save_smite never removes a folder, and replaces a file as a write would" begin
    cam = IdealCamera(1:64, 1:64, 0.1)
    mk(x) = SmiteSMLD{Float64,Emitter2DFit{Float64}}(
        [Emitter2DFit{Float64}(x, 1.5, 1000.0, 10.0, 0.01, 0.01, 50.0, 2.0)], cam, 1, 1, Dict{String,Any}())
    xof(dir, name) = load_smite_2d(SmiteSMD(dir, name)).emitters[1].x
    mktempdir() do dir
        # A target that is a folder, or a name that is no file: an error, and every file stays
        mkdir(joinpath(dir, "Cell_01")); write(joinpath(dir, "Cell_01", "keep.txt"), "keep")
        @test_throws r"is a folder" save_smite(mk(1.0), dir, "Cell_01")
        for name in ("", ".", "..")
            @test_throws r"names no file" save_smite(mk(1.0), dir, name)
        end
        @test read(joinpath(dir, "Cell_01", "keep.txt"), String) == "keep"
        @test readdir(dir) == ["Cell_01"]

        # A missing folder, or a folder argument that is a file, is named
        err = try save_smite(mk(1.0), joinpath(dir, "missing"), "r.mat"); nothing catch ex; ex end
        @test err isa ArgumentError && occursin("missing", err.msg) && occursin("does not exist", err.msg)
        write(joinpath(dir, "afile"), "x")
        @test_throws r"is not a folder" save_smite(mk(1.0), joinpath(dir, "afile"), "r.mat")
        rm(joinpath(dir, "afile"))

        # A file name with a folder in it, and an absolute file name, save as before
        save_smite(mk(2.0), dir, joinpath("Cell_01", "r.mat"))
        @test xof(joinpath(dir, "Cell_01"), "r.mat") ≈ 2.0
        save_smite(mk(3.0), dir, joinpath(dir, "abs.mat"))
        @test xof(dir, "abs.mat") ≈ 3.0

        # An existing file is rewritten in place, as a direct write would: it keeps its
        # permissions and its hard links, and a file the system will not let us write is refused
        # and unchanged (POSIX permissions; CI runs on Linux)
        if !Sys.iswindows()
            f = joinpath(dir, "mode.mat")
            save_smite(mk(1.0), dir, "mode.mat"); chmod(f, 0o640)
            hardlink(f, joinpath(dir, "hard.mat"))
            save_smite(mk(4.0), dir, "mode.mat")
            @test filemode(f) & 0o777 == 0o640 && xof(dir, "mode.mat") ≈ 4.0
            @test xof(dir, "hard.mat") ≈ 4.0   # the same file, not a replaced one
            chmod(f, 0o440)
            if Libc.geteuid() != 0   # root may write any file
                @test_throws r"cannot be written" save_smite(mk(5.0), dir, "mode.mat")
                @test xof(dir, "mode.mat") ≈ 4.0
            end
            chmod(f, 0o640)

            # A symlink is followed: the linked file is rewritten and the link stays a link; a
            # symlink to a folder is a folder
            save_smite(mk(1.0), dir, "real.mat")
            symlink(joinpath(dir, "real.mat"), joinpath(dir, "link.mat"))
            save_smite(mk(6.0), dir, "link.mat")
            @test islink(joinpath(dir, "link.mat"))
            @test xof(dir, "real.mat") ≈ 6.0
            symlink(joinpath(dir, "Cell_01"), joinpath(dir, "dirlink"))
            @test_throws r"is a folder" save_smite(mk(1.0), dir, "dirlink")
            @test read(joinpath(dir, "Cell_01", "keep.txt"), String) == "keep"

            # A symlink to nothing is refused, and the link stays
            symlink(joinpath(dir, "made.mat"), joinpath(dir, "dangling.mat"))
            @test_throws r"symlink to nothing" save_smite(mk(7.0), dir, "dangling.mat")
            @test islink(joinpath(dir, "dangling.mat")) && !ispath(joinpath(dir, "made.mat"))

            # Something that is not a file is refused and stays
            run(`mkfifo $(joinpath(dir, "fifo"))`)
            @test_throws r"not a regular file" save_smite(mk(1.0), dir, "fifo")
            @test ispath(joinpath(dir, "fifo"))
        end
        # Nothing else is left in the folder
        @test issubset(readdir(dir), ["Cell_01", "abs.mat", "mode.mat", "hard.mat", "real.mat",
                                      "link.mat", "dirlink", "dangling.mat", "fifo"])
    end
end

@testset "paths of any string type; SmiteSMD shows its whole path" begin
    cam = IdealCamera(1:64, 1:64, 0.1)
    e = Emitter2DFit{Float64}(1.0, 1.5, 1000.0, 10.0, 0.01, 0.01, 50.0, 2.0)
    smld = SmiteSMLD{Float64,Emitter2DFit{Float64}}([e], cam, 1, 1, Dict{String,Any}())
    mktempdir() do dir
        sub(s) = SubString(string("xx", s, "yy"), 3, 2 + length(s))
        save_smite(smld, sub(dir), sub("sub.mat"))
        @test isfile(joinpath(dir, "sub.mat"))
        smd = SmiteSMD(sub(dir), sub("sub.mat"))
        @test smd isa SmiteSMD && smd.filepath isa String && smd.filename isa String
        @test SmiteSMD(sub(dir), sub("sub.mat"), sub("SMD")).varname == "SMD"
        @test load_smite_2d(smd).emitters[1].x ≈ 1.0
        # Mutable as before
        smd.filename = sub("other.mat")
        @test smd.filename == "other.mat" && smd.filename isa String
    end

    smd = SmiteSMD("/a/b/c", "x.mat")
    @test repr(smd) == "SmiteSMD(\"/a/b/c\", \"x.mat\", \"SMD\")"
    @test repr(eval(Meta.parse(repr(smd)))) == repr(smd)   # it reads as the constructor call it is
    @test occursin("/a/b/c", sprint(show, MIME("text/plain"), smd))
end

@testset "save_smite applies the loaders' rule to counts and indices" begin
    cam = IdealCamera(1:64, 1:64, 0.1)
    two53 = Int64(2)^53
    big = two53 + 1                   # the first whole number a double cannot hold
    sm(es, nf, nd) = AnyCountSMLD(es, cam, nf, nd, Dict{String,Any}())
    a = _anyloc()
    mktempdir() do dir
        # What it writes loads: n_frames 0 ("unknown") and n_datasets 0 included
        for (nf, nd) in ((0, 0), (0, 3), (5, 0))
            save_smite(SMLMCore.BasicSMLD([Emitter2DFit{Float64}(1.0, 1.5, 1000.0, 10.0, 0.01, 0.01, 50.0, 2.0)],
                                          cam, nf, nd), dir, "zero.mat")
            back = load_smite_2d(SmiteSMD(dir, "zero.mat"))
            @test back.n_frames == nf && back.n_datasets == nd
        end
        # Whole numbers a double holds exactly, of any number type, write and load back exactly
        # (2^53 + 2 is such a double)
        if Sys.WORD_SIZE == 64
            save_smite(sm([_anyloc(frame=two53, dataset=UInt64(5), track_id=two53 + 2)], two53 + 2, 7 // 1),
                       dir, "edge.mat")
            back = load_smite_2d(SmiteSMD(dir, "edge.mat"))
            @test back.n_frames == two53 + 2 && back.n_datasets == 7
            e = back.emitters[1]
            @test e.frame == two53 && e.dataset == 5 && e.track_id == two53 + 2
        end
        # Anything the loaders would refuse or read back differently is an ArgumentError naming
        # the value; nothing is written
        for (name, smld, msg) in (
                ("nf.mat", sm([a], big, 1), "n_frames is $big"),
                ("nfr.mat", sm([a], (big)//1, 1), "n_frames is $(repr(big//1))"),
                ("nfu.mat", sm([a], UInt64(1) << 63, 1), "n_frames is $(repr(UInt64(1) << 63))"),
                ("nfneg.mat", sm([a], -1, 1), "n_frames is -1; SMITE stores it as a MATLAB double, so it must be a whole number of at least 0"),
                ("nfh.mat", sm([a], 2.5, 1), "n_frames is 2.5"),
                ("nd.mat", sm([a], 1, big), "n_datasets is $big"),
                ("ndneg.mat", sm([a], 1, -2), "n_datasets is -2"),
                ("fr.mat", sm([a, _anyloc(frame=big)], 1, 1), "frame of emitter 2 is $big"),
                ("frh.mat", sm([a, a, _anyloc(frame=2.5)], 1, 1), "frame of emitter 3 is 2.5"),
                ("ds.mat", sm([_anyloc(dataset=big)], 1, 1), "dataset of emitter 1 is $big"),
                ("dsu.mat", sm([_anyloc(dataset=typemax(UInt64))], 1, 1), "dataset of emitter 1 is $(repr(typemax(UInt64)))"),
                ("tr.mat", sm([a, a, _anyloc(track_id=-big)], 1, 1), "track_id of emitter 3 is -$big"),
                ("trm.mat", sm([_anyloc(track_id=missing)], 1, 1), "track_id of emitter 1 is missing"))
            err = try save_smite(smld, dir, name); nothing catch ex; ex end
            @test err isa ArgumentError
            err isa ArgumentError && @test occursin(msg, err.msg) && occursin("MATLAB double", err.msg)
            @test !isfile(joinpath(dir, name))
        end
        # An existing file is left as it was
        save_smite(sm([a], 1, 1), dir, "keep.mat")
        before = read(joinpath(dir, "keep.mat"))
        @test_throws ArgumentError save_smite(sm([a], big, 1), dir, "keep.mat")
        @test read(joinpath(dir, "keep.mat")) == before
    end
end

@testset "save_smite writes every emitter column exactly, as single or double" begin
    cam = IdealCamera(1:64, 1:64, 0.1)
    sm(es) = AnyCountSMLD(es, cam, 1, 1, Dict{String,Any}())
    mktempdir() do dir
        # Integers and Float16 values become doubles exactly (MAT cannot write Float16)
        save_smite(sm([_anyloc(photons=Int64(1000), bg=Float16(2.5)), _anyloc(photons=1200.0, bg=Float16(3))]),
                   dir, "ints.mat")
        smd = matread(joinpath(dir, "ints.mat"))["SMD"]
        @test smd["Photons"] isa AbstractArray{Float64} && vec(smd["Photons"]) == [1000.0, 1200.0]
        @test smd["Bg"] isa AbstractArray{Float64} && vec(smd["Bg"]) == [2.5, 3.0]
        @test [e.photons for e in load_smite_2d(SmiteSMD(dir, "ints.mat")).emitters] == [1000.0, 1200.0]
        # A complex value is written as a complex double column, as 0.1.0 wrote it; the loaders
        # drop that localization, as they do for SMITE's complex fits
        save_smite(sm([_anyloc(photons=1.0 + 1.0im), _anyloc(photons=900.0)]), dir, "pc.mat")
        @test matread(joinpath(dir, "pc.mat"))["SMD"]["Photons"] isa AbstractArray{ComplexF64}
        back = @test_logs (:warn, r"non-zero imaginary") load_smite_2d(SmiteSMD(dir, "pc.mat"))
        @test [e.photons for e in back.emitters] == [900.0]
        # A number a double would round, or a value that is not a number, is an ArgumentError
        # naming the emitter (0.1.0 wrote a cell its loader refused); nothing is written
        for (name, es, msg) in (
                ("p53.mat", [_anyloc(photons=Int64(2)^53 + 1), _anyloc(photons=1.0)], "photons of emitter 1 is 9007199254740993"),
                ("bgbig.mat", [_anyloc(), _anyloc(bg=big"0.1")], "bg of emitter 2 is"),
                ("pmiss.mat", [_anyloc(photons=missing)], "photons of emitter 1 is missing"),
                ("xbig.mat", [_anyloc(x=big"0.1")], "x of emitter 1 is"),
                ("pcbig.mat", [_anyloc(photons=big"0.1" + 1im)], "photons of emitter 1 is"),
                ("ptext.mat", [_anyloc(), _anyloc(bg="10")], "bg of emitter 2 is \"10\""))
            err = try save_smite(sm(es), dir, name); nothing catch ex; ex end
            @test err isa ArgumentError
            err isa ArgumentError && @test occursin(msg, err.msg) && occursin("MATLAB double holds exactly", err.msg)
            @test !isfile(joinpath(dir, name))
        end
    end
end

@testset "save_smite checks the camera it writes and computes positions in Float64" begin
    em = Emitter2DFit{Float64}(1.0, 1.5, 1000.0, 10.0, 0.01, 0.01, 50.0, 2.0)
    mktempdir() do dir
        # A camera the loaders would refuse is an ArgumentError, before anything is written
        for (name, cam, msg) in (
                ("ny0.mat", IdealCamera(2, 0, 0.125), "the camera has 2 × 0 pixels"),
                ("wide.mat", IdealCamera(2^24 + 1, 1, 1.0f0), "16777217 × 1 pixels of 1.0 µm do not give finite, increasing"),
                ("inf.mat", IdealCamera(1, 1, big"1e309"), "the camera's pixel size Inf is not finite"))
            err = try save_smite(SMLMCore.BasicSMLD([em], cam, 1, 1), dir, name); nothing catch ex; ex end
            @test err isa ArgumentError
            err isa ArgumentError && @test occursin(msg, err.msg)
            @test !isfile(joinpath(dir, name))
        end
        # Positions in SMITE's pixels are computed in Float64 from the written pixel size: a
        # position a double holds loads back exactly, with a Float32 camera too
        cam32 = IdealCamera(1:64, 1:64, 0.125f0)
        x = 33554433 // 33554432                       # 1 + 2^-25 µm
        save_smite(AnyCountSMLD([_anyloc(x=x)], cam32, 1, 1, Dict{String,Any}()), dir, "exact.mat")
        @test matread(joinpath(dir, "exact.mat"))["SMD"]["X"] === 8.500000238418579
        @test load_smite_2d(SmiteSMD(dir, "exact.mat")).emitters[1].x === Float64(x)
        # A finite position whose pixel coordinate a single cannot hold is an ArgumentError
        e32 = Emitter2DFit{Float32}(floatmax(Float32), 1.0f0, 1000.0f0, 10.0f0, 0.01f0, 0.01f0, 50.0f0, 2.0f0)
        @test_throws r"x of emitter 1 is .*a MATLAB single in SMITE's pixels cannot hold" save_smite(
            SMLMCore.BasicSMLD([e32], IdealCamera(1:64, 1:64, 0.25f0), 1, 1), dir, "over.mat")
        @test !isfile(joinpath(dir, "over.mat"))
        # A finite value the loaders would read back as infinite, or an uncertainty that would
        # read back as 0, is an ArgumentError too
        cam3 = IdealCamera(1:64, 1:64, 3.0)
        for (name, es, msg) in (("xmax.mat", [_anyloc(x=floatmax(Float64))], "x of emitter 1 is"),
                                ("sxmax.mat", [AnyNumLoc(1.0, 1.5, 1000.0, 10.0, floatmax(Float64), 0.01, 50.0, 2.0, 1, 1, 0, 1)],
                                 "σ_x of emitter 1 is"))
            @test_throws Regex(msg) save_smite(AnyCountSMLD(es, cam3, 1, 1, Dict{String,Any}()), dir, name)
            @test !isfile(joinpath(dir, name))
        end
        tiny = Emitter2DFit{Float32}(1.0f0, 1.0f0, 1000.0f0, 10.0f0, nextfloat(0.0f0), 0.01f0, 50.0f0, 2.0f0)
        @test_throws r"σ_x of emitter 1 is" save_smite(
            SMLMCore.BasicSMLD([tiny], IdealCamera(1:64, 1:64, 2.0f0), 1, 1), dir, "tiny.mat")
        # A camera origin that is not a double is subtracted exactly
        big0 = BigFloat(2)^53 + 1
        camb = IdealCamera(collect(big0 .+ (0:4)), collect(BigFloat(0):4))
        @test_logs (:warn, r"first pixel edge is at") save_smite(
            AnyCountSMLD([_anyloc(x=Float64(2)^53 + 2)], camb, 1, 1, Dict{String,Any}()), dir, "bigorigin.mat")
        @test matread(joinpath(dir, "bigorigin.mat"))["SMD"]["X"] === 1.5
        # A camera whose edges are not those the loaders rebuild from the written pixel size is
        # refused: equal edges, a Float32 camera whose stored spacing far from the origin is not
        # its first pixel's, a Float16 camera (its edges round differently from the double written)
        for (name, cam) in (("flat.mat", IdealCamera(Float32[0, 1], Float32[2^24, 2^24])),
                            ("far32.mat", IdealCamera(Float32[0, 0.1], Float32[10000, 10000.1])),
                            ("f16.mat", IdealCamera(512, 1, Float16(0.1))))
            @test_throws r"not uniform and square" save_smite(SMLMCore.BasicSMLD([em], cam, 1, 1), dir, name)
            @test !isfile(joinpath(dir, name))
        end
        # Edges are subtracted exactly: Float16 edges offset from a pixel boundary round onto
        # another grid (refused), and an exact Int8 grid away from 0 does not overflow (saved)
        p16 = Float16(0.125)
        @test_throws r"not uniform and square" save_smite(SMLMCore.BasicSMLD([em],
            IdealCamera(Float16[0, p16], Float16.(0:65) .* p16 .- Float16(0.003)), 1, 1), dir, "f16off.mat")
        int8 = IdealCamera(collect(Int8(-100):Int8(100)), Int8[0, 1])
        e8 = Emitter2DFit{Float64}(-50.0, 0.5, 1000.0, 10.0, 0.01, 0.01, 50.0, 2.0)
        @test_logs (:warn, r"first pixel edge is at") save_smite(SMLMCore.BasicSMLD([e8], int8, 1, 1), dir, "int8.mat")
        smd = matread(joinpath(dir, "int8.mat"))["SMD"]
        @test smd["XSize"] == 200 && smd["YSize"] == 1 && smd["PixelSize"] === 1.0 && smd["X"] === 50.5
        # A large Float32 camera at 0 is rebuilt exactly, and saves
        save_smite(SMLMCore.BasicSMLD([em], IdealCamera(4096, 4096, 0.1f0), 1, 1), dir, "big32.mat")
        @test load_smite_2d(SmiteSMD(dir, "big32.mat")).camera.pixel_edges_x == IdealCamera(4096, 4096, 0.1f0).pixel_edges_x
        # Values with a zero imaginary part are their real part, as 0.1.0 saved them
        save_smite(AnyCountSMLD([_anyloc(photons=4.0 + 0.0im)], IdealCamera(1:64, 1:64, 0.1), 1, 1,
                                Dict{String,Any}()), dir, "c0.mat")
        @test load_smite_2d(SmiteSMD(dir, "c0.mat")).emitters[1].photons === 4.0
    end
end

@testset "save_smite writes Z one way for emitters and metadata" begin
    cam = IdealCamera(1:64, 1:64, 0.1)
    sm(v; md=Dict{String,Any}()) = SmiteSMLD{Float64,eltype(v)}(v, cam, 1, 1, md)
    A = SMLMCore.AbstractEmitter
    e2 = Emitter2DFit{Float64}(1.0, 1.5, 1000.0, 10.0, 0.01, 0.01, 50.0, 2.0)
    e3 = Emitter3DFit{Float64}(1.0, 1.5, 0.75, 1000.0, 10.0, 0.01, 0.01, 0.125, 50.0, 2.0)
    mktempdir() do dir
        # Complex z is written as a complex column, as 0.1.0 did; load_smite_3d drops that row
        save_smite(sm([_oddz(0.25 + 1.0im, 0.02), _oddz(0.5 + 0.0im, 0.02)]), dir, "cz.mat")
        smd = matread(joinpath(dir, "cz.mat"))["SMD"]
        @test smd["NDims"] == 3 && smd["Z"] isa AbstractArray{ComplexF64}
        back = @test_logs (:warn, r"non-zero imaginary") load_smite_3d(SmiteSMD(dir, "cz.mat"))
        @test [e.z for e in back.emitters] == [0.5]
        # A number a double would round is an error wherever it is in the column, also after a
        # value that is not a number (which alone saves without Z)
        for zs in ([missing, big"0.1"], [big"0.1", missing])
            @test_throws r"z of emitter \d is 0.1" save_smite(sm([_oddz(z, 0.02) for z in zs]), dir, "zo.mat")
            @test_throws r"σ_z of emitter \d is 0.1" save_smite(sm([_oddz(0.5, z) for z in zs]), dir, "zo.mat")
        end
        @test !isfile(joinpath(dir, "zo.mat"))
        # Metadata Z and Z_SE never stand in for emitters whose own Z was not saved
        md() = Dict{String,Any}("Z" => [10.0, 20.0], "Z_SE" => [0.1, 0.1])
        @test_logs (:warn, r"mixed 2D/3D") save_smite(sm(A[e2, e3]; md=md()), dir, "mixz.mat")
        @test_logs (:warn, r"3D emitters") save_smite(sm([_oddz(missing, 0.02), _oddz(missing, 0.02)]; md=md()),
                                                     dir, "missz.mat")
        for f in ("mixz.mat", "missz.mat")
            smd = matread(joinpath(dir, f))["SMD"]
            @test !haskey(smd, "Z") && !haskey(smd, "Z_SE") && smd["NDims"] == 2
        end
        # On 2D emitters, metadata Z and Z_SE with one value per emitter (an SMLD built by hand, or
        # one SMLMData 0.7.1 loaded from a 3D file) are written as the emitters' z would be, and a
        # save of the saved file's load agrees with the first save
        twice(name, mdict) = begin
            save_smite(sm([e2, e2]; md=mdict), dir, name)
            once = matread(joinpath(dir, name))["SMD"]
            # (the load warns that it leaves out Z, which the file now holds for each localization)
            back = Base.CoreLogging.with_logger(() -> load_smite_2d(SmiteSMD(dir, name)),
                                                Base.CoreLogging.NullLogger())
            save_smite(back, dir, "again-" * name)
            once, matread(joinpath(dir, "again-" * name))["SMD"]
        end
        for (a, b) in (twice("mdz.mat", Dict{String,Any}("Z" => [0.25, 0.5], "Z_SE" => [0.1, 0.2])),
                       twice("mdz3.mat", Dict{String,Any}("Z" => [0.25, 0.5], "Z_SE" => [0.1, 0.2], "NDims" => 3.0)))
            @test a["NDims"] == 3 && b["NDims"] == 3
            @test vec(a["Z"]) == vec(b["Z"]) == [0.25, 0.5]
        end
        # (the first save drops them with a warning, so the saved file has none to drop again)
        a, b = @test_logs (:warn, r"Z, Z_SE not saved") twice(
            "mdzbad.mat", Dict{String,Any}("Z" => Any["bad", 2.0], "Z_SE" => [0.1, 0.2]))
        @test !haskey(a, "Z") && !haskey(b, "Z") && a["NDims"] == b["NDims"] == 2
        @test_throws r"Z of metadata row 1 is 9007199254740993" save_smite(
            sm([e2, e2]; md=Dict{String,Any}("Z" => Int64[2^53 + 1, 1], "Z_SE" => [0.1, 0.2])), dir, "mdzbig.mat")
        # Data whose emitters' ids are rows of the metadata's Z (SMLMData 0.7.1 set id to the file
        # row and kept a 3D file's Z in file order): each emitter gets its own row's z, after a
        # reordering or a localization dropped by the complex screening alike
        f(fr, id) = Emitter2DFit{Float64}(1.0, 1.5, 1000.0, 10.0, 0.01, 0.01, 50.0, 2.0; frame=fr, id=id)
        zmd() = Dict{String,Any}("Z" => [0.25, 0.5, 0.75], "Z_SE" => [0.1, 0.2, 0.3])
        for (es, z, σz) in (([f(3, 3), f(1, 1)], [0.75, 0.25], [0.3, 0.1]),
                            ([f(1, 1), f(3, 3)], [0.25, 0.75], [0.1, 0.3]),
                            ([f(2, 2), f(1, 1), f(3, 3)], [0.5, 0.25, 0.75], [0.2, 0.1, 0.3]))
            save_smite(sm(es; md=zmd()), dir, "byid.mat")
            back = load_smite_3d(SmiteSMD(dir, "byid.mat"))
            @test [e.frame for e in back.emitters] == [e.frame for e in es]
            @test [e.z for e in back.emitters] == z && [e.σ_z for e in back.emitters] ≈ σz
        end
        # Ids that are not rows of a Z with as many values as emitters: not taken for coordinates,
        # dropped with a warning (a load would take them for each emitter's)
        @test_logs (:warn, r"Z, Z_SE not saved: their values cannot be matched") save_smite(
            sm([f(1, 5), f(2, 6)]; md=Dict{String,Any}("Z" => [0.25, 0.5], "Z_SE" => [0.1, 0.2])), dir, "noid.mat")
        smd = matread(joinpath(dir, "noid.mat"))["SMD"]
        @test !haskey(smd, "Z") && smd["NDims"] == 2
        # A tuple in the metadata is the array it holds
        save_smite(sm([e2, e2]; md=Dict{String,Any}("Z" => (0.25, 0.5), "Z_SE" => [0.1, 0.2])), dir, "tup.mat")
        smd = matread(joinpath(dir, "tup.mat"))["SMD"]
        @test smd["NDims"] == 3 && vec(smd["Z"]) == [0.25, 0.5]
        @test_throws r"Z of metadata row 1 is 9007199254740993" save_smite(
            sm([e2, e2]; md=Dict{String,Any}("Z" => (Int64(2)^53 + 1, 1), "Z_SE" => [0.1, 0.2])), dir, "tupbig.mat")
        # Metadata Z of another length is ordinary metadata, written as given; the file stays 2D
        save_smite(sm([e2, e2]; md=Dict{String,Any}("Z" => [1.0, 2.0, 3.0])), dir, "mdz3rows.mat")
        smd = matread(joinpath(dir, "mdz3rows.mat"))["SMD"]
        @test smd["NDims"] == 2 && vec(smd["Z"]) == [1.0, 2.0, 3.0]
    end
end
