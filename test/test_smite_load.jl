# SMITE loading. Each fixture has the layout of a SMITE results file: per-localization columns
# n x 1, scalars Float64, fields SMITE did not compute 0 x 0.
using SMLMio, Test
using SMLMCore: SMLMCore, IdealCamera, BasicSMLD, Emitter2DFit, Emitter3DFit, get_pixel_centers
using MAT: matread, matwrite

col(v) = reshape(collect(v), :, 1)

function smite_smd(n::Int=3; pixel_size=0.1, xsize=64.0, ysize=64.0)
    return Dict{String,Any}(
        "X" => col(Float32.(range(10, 40; length=n))), "Y" => col(Float32.(range(5, 50; length=n))),
        "X_SE" => col(fill(0.2f0, n)), "Y_SE" => col(fill(0.3f0, n)),
        "Photons" => col(fill(1000.0f0, n)), "Bg" => col(fill(5.0f0, n)),
        "Photons_SE" => col(fill(30.0f0, n)), "Bg_SE" => col(fill(0.5f0, n)),
        "FrameNum" => col(Float64.(1:n)), "DatasetNum" => col(ones(n)), "ConnectID" => col(Float64.(1:n)),
        "PValue" => col(fill(0.5f0, n)), "ThreshFlag" => col(zeros(UInt32, n)),
        "DriftX" => zeros(Float32, 100, 1), "Z" => zeros(0, 0), "Z_SE" => zeros(0, 0), "ZOffset" => zeros(0, 0),
        "NDims" => 2.0, "NFrames" => 100.0, "NDatasets" => 1.0, "PixelSize" => pixel_size,
        "XSize" => xsize, "YSize" => ysize)
end

# The same with SMITE's 3D columns: Z and Z_SE in microns, NDims 3.
function smite_smd_3d(n::Int=3; kwargs...)
    s = smite_smd(n; kwargs...)
    s["Z"] = col(Float32.(range(-0.4, 0.4; length=n)))
    s["Z_SE"] = col(fill(0.04f0, n))
    s["NDims"] = 3.0
    return s
end

write_smd(dir, name, s) = (matwrite(joinpath(dir, name), Dict("SMD" => s)); SmiteSMD(dir, name))
quiet(f) = Base.CoreLogging.with_logger(f, Base.CoreLogging.NullLogger())

@testset "SMITE loading" begin
    @testset "pixels load as microns" begin
        # SMITE: X, Y, X_SE, Y_SE in pixels, (1, 1) the center of the top-left pixel; Z, Z_SE microns.
        mktempdir() do dir
            px = 0.108
            for three_d in (false, true)
                s = three_d ? smite_smd_3d(3; pixel_size=px) : smite_smd(3; pixel_size=px)
                three_d && (s["ZSize"] = 1.0)   # present here; its absence is tested below
                smd = write_smd(dir, "units$three_d.mat", s)
                smld = three_d ? load_smite_3d(smd) : load_smite_2d(smd)
                es = smld.emitters
                @test [e.x for e in es] ≈ (vec(s["X"]) .- 0.5) .* px
                @test [e.y for e in es] ≈ (vec(s["Y"]) .- 0.5) .* px
                @test [e.σ_x for e in es] ≈ vec(s["X_SE"]) .* px
                @test [e.σ_y for e in es] ≈ vec(s["Y_SE"]) .* px
                @test [e.photons for e in es] ≈ vec(s["Photons"])
                if three_d
                    @test [e.z for e in es] ≈ vec(s["Z"])
                    @test [e.σ_z for e in es] ≈ vec(s["Z_SE"])
                end
                edges = smld.camera.pixel_edges_x
                @test all(e -> edges[1] <= e.x <= edges[end], es)
            end
            # The center of SMITE's pixel (1, 1) is the center of the camera's first pixel.
            s = smite_smd(2; pixel_size=px)
            s["X"] = col(Float32[1, 64]); s["Y"] = col(Float32[1, 64])
            smld = load_smite_2d(write_smd(dir, "corner.mat", s))
            cx, cy = get_pixel_centers(smld.camera)
            @test smld.emitters[1].x ≈ cx[1]
            @test smld.emitters[1].y ≈ cy[1]
            @test smld.emitters[2].x ≈ cx[end]
            @test smld.emitters[2].y ≈ cy[end]
        end
    end

    @testset "save_smite writes SMITE pixels" begin
        mktempdir() do dir
            # Microns in, pixels out: the center of the camera's first pixel is SMITE's (1, 1).
            px = 0.108
            cam = IdealCamera(1:64, 1:64, px)
            es = [Emitter2DFit{Float64}(0.5px, 0.5px, 1000.0, 5.0, 0.02, 0.03, 30.0, 0.5),
                  Emitter2DFit{Float64}(3.0, 4.0, 1000.0, 5.0, 0.02, 0.03, 30.0, 0.5)]
            smld = SmiteSMLD{Float64,Emitter2DFit{Float64}}(es, cam, 1, 1, Dict{String,Any}())
            save_smite(smld, dir, "out.mat")
            out = matread(joinpath(dir, "out.mat"))["SMD"]
            @test vec(out["X"]) ≈ [1.0, 3.0 / px + 0.5]
            @test vec(out["Y"]) ≈ [1.0, 4.0 / px + 0.5]
            @test vec(out["X_SE"]) ≈ [0.02, 0.02] ./ px
            @test vec(out["Y_SE"]) ≈ [0.03, 0.03] ./ px
            # A SMITE file saved after loading holds its original pixel values.
            s = smite_smd(3; pixel_size=px)
            loaded = load_smite_2d(write_smd(dir, "in.mat", s))
            save_smite(loaded, dir, "again.mat")
            again = matread(joinpath(dir, "again.mat"))["SMD"]
            for k in ("X", "Y", "X_SE", "Y_SE")
                @test vec(again[k]) ≈ vec(s[k]) rtol = 1e-6
            end
        end
    end
end
