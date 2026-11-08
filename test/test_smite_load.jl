# SMITE loading. Each fixture has the layout of a SMITE results file: per-localization columns
# n x 1, scalars Float64, fields SMITE did not compute 0 x 0.
using SMLMio, Test
using SMLMCore: SMLMCore, IdealCamera, BasicSMLD, Emitter2DFit, Emitter3DFit, get_pixel_centers
using MAT: matread, matwrite

col(v) = reshape(collect(v), :, 1)

function smite_smd(n::Int=3; pixel_size=0.1, xsize=64.0, ysize=64.0)
    return Dict{String,Any}(
        "X" => col(Float32.(10 .+ 30 .* (0:n-1) ./ max(n - 1, 1))), "Y" => col(Float32.(5 .+ 45 .* (0:n-1) ./ max(n - 1, 1))),
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
    s["Z"] = col(Float32.(-0.4 .+ 0.8 .* (0:n-1) ./ max(n - 1, 1)))
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

    @testset "fields without one value per localization" begin
        mktempdir() do dir
            # Older SMITE versions wrote Photons_SE and Bg_SE with one value per dataset:
            # they load as NaN with a warning, raw values kept in metadata.
            s = smite_smd(5)
            s["NDatasets"] = 2.0
            s["Photons_SE"] = col(Float32[31, 32]); s["Bg_SE"] = col(Float32[0.6, 0.7])
            smd = write_smd(dir, "perdataset.mat", s)
            @test_logs (:warn, r"Photons_SE, Bg_SE do not hold one value per localization") load_smite_2d(smd)
            smld = quiet(() -> load_smite_2d(smd))
            @test length(smld.emitters) == 5
            @test all(e -> isnan(e.σ_photons) && isnan(e.σ_bg), smld.emitters)
            @test vec(smld.metadata["Photons_SE"]) == Float32[31, 32]
            @test vec(smld.metadata["Bg_SE"]) == Float32[0.6, 0.7]
            @test smld.metadata["smite_defaulted_fields"] == ["Photons_SE", "Bg_SE"]
            @test [e.photons for e in smld.emitters] ≈ vec(s["Photons"])
            # Empty (not computed) Photons_SE and Bg_SE: NaN, no warning.
            s = smite_smd(4)
            s["Photons_SE"] = zeros(0, 0); s["Bg_SE"] = zeros(0, 0)
            smld = @test_logs load_smite_2d(write_smd(dir, "emptyse.mat", s))
            @test all(e -> isnan(e.σ_photons) && isnan(e.σ_bg), smld.emitters)
            # ConnectID is empty before frame connection: track_id 0, no warning.
            s = smite_smd(4)
            s["ConnectID"] = zeros(UInt32, 0, 0)
            smld = @test_logs load_smite_2d(write_smd(dir, "noconnect.mat", s))
            @test [e.track_id for e in smld.emitters] == zeros(Int, 4)
            @test smld.metadata["smite_defaulted_fields"] == ["ConnectID"]
            # No DatasetNum: dataset 1.
            s = smite_smd(3)
            delete!(s, "DatasetNum")
            smld = load_smite_2d(write_smd(dir, "nodataset.mat", s))
            @test [e.dataset for e in smld.emitters] == [1, 1, 1]
            # The same in 3D.
            s = smite_smd_3d(3)
            s["ZSize"] = 1.0
            s["ConnectID"] = zeros(UInt32, 0, 0); s["Bg_SE"] = zeros(0, 0)
            smld = load_smite_3d(write_smd(dir, "noconnect3d.mat", s))
            @test [e.track_id for e in smld.emitters] == [0, 0, 0]
            @test all(e -> isnan(e.σ_bg), smld.emitters)
            # A required field of the wrong length is an error naming it.
            s = smite_smd(3)
            s["Photons"] = col(Float32[1000, 1000])
            @test_throws r"Photons holds 2 values for 3 localizations" load_smite_2d(write_smd(dir, "short.mat", s))
        end
    end

    @testset "pixel size from the file or the caller, never assumed" begin
        mktempdir() do dir
            # No PixelSize, or an empty one: the caller passes pixel_size; nothing is assumed.
            for (name, pxfield) in (("absent", nothing), ("empty", zeros(0, 0)))
                s = smite_smd(3)
                pxfield === nothing ? delete!(s, "PixelSize") : (s["PixelSize"] = pxfield)
                smd = write_smd(dir, "px$name.mat", s)
                @test_throws r"no PixelSize; pass pixel_size" load_smite_2d(smd)
                smld = load_smite_2d(smd; pixel_size=0.108)
                @test [e.x for e in smld.emitters] ≈ (vec(s["X"]) .- 0.5) .* 0.108
                @test smld.camera.pixel_edges_x[end] ≈ 64 * 0.108
                @test smld.metadata["pixel_size"] == 0.108
                s3 = smite_smd_3d(3); s3["ZSize"] = 1.0
                pxfield === nothing ? delete!(s3, "PixelSize") : (s3["PixelSize"] = pxfield)
                smld = load_smite_3d(write_smd(dir, "px3$name.mat", s3); pixel_size=0.108)
                @test [e.y for e in smld.emitters] ≈ (vec(s3["Y"]) .- 0.5) .* 0.108
            end
            # The file's PixelSize wins; a different pixel_size is an error, the same one is fine.
            smd = write_smd(dir, "pxfile.mat", smite_smd(3; pixel_size=0.0978))
            @test load_smite_2d(smd; pixel_size=0.0978).camera.pixel_edges_x[2] ≈ 0.0978
            @test_throws r"PixelSize 0.0978 µm, but pixel_size = 0.1" load_smite_2d(smd; pixel_size=0.1)
            s = smite_smd(3); delete!(s, "PixelSize")
            @test_throws r"not positive" load_smite_2d(write_smd(dir, "pxneg.mat", s); pixel_size=-0.1)
        end
    end

    @testset "3D without ZSize" begin
        # SMITE's createSMD defines XSize and YSize but no ZSize. This SMD is hand-built with
        # SMITE's fields.
        mktempdir() do dir
            s = smite_smd_3d(3)
            @test !haskey(s, "ZSize")
            smld = load_smite_3d(write_smd(dir, "nozsize.mat", s))
            @test length(smld.emitters) == 3
            @test [e.z for e in smld.emitters] ≈ vec(s["Z"])
            @test smld.metadata["data_size"] == [64, 64]
            s["ZSize"] = 7.0
            @test load_smite_3d(write_smd(dir, "zsize.mat", s)).metadata["data_size"] == [64, 64, 7]
        end
    end

    @testset "results without localizations" begin
        mktempdir() do dir
            # SMITE leaves XSize, YSize, NFrames, NDatasets and every column empty in a result
            # with no localizations: no camera can be built, and the error says so.
            s = smite_smd(1)
            for (k, v) in s
                (v isa Number || v isa AbstractArray) && k != "PixelSize" && (s[k] = zeros(0, 0))
            end
            @test_throws r"no image size \(XSize, YSize absent or empty; 0 localizations\)" load_smite_2d(write_smd(dir, "empty.mat", s))
            # No localizations but an image size: an empty data set with its camera.
            s = smite_smd(1)
            for k in ("X", "Y", "X_SE", "Y_SE", "Photons", "Bg", "Photons_SE", "Bg_SE", "FrameNum",
                      "DatasetNum", "ConnectID", "PValue", "ThreshFlag")
                s[k] = zeros(0, 0)
            end
            smld = load_smite_2d(write_smd(dir, "none.mat", s))
            @test isempty(smld.emitters)
            @test length(smld.camera.pixel_edges_x) == 65
            # No NFrames: an error that names it.
            s = smite_smd(2); delete!(s, "NFrames")
            @test_throws r"no NFrames or NDatasets" load_smite_2d(write_smd(dir, "nonframes.mat", s))
        end
    end

    @testset "row vectors and a single localization" begin
        mktempdir() do dir
            # Real files store columns; an SMD assembled in MATLAB may hold rows.
            s = smite_smd(4)
            for (k, v) in s
                v isa AbstractMatrix && size(v, 2) == 1 && size(v, 1) == 4 && (s[k] = permutedims(v))
            end
            smld = load_smite_2d(write_smd(dir, "rows.mat", s))
            @test length(smld.emitters) == 4
            @test [e.frame for e in smld.emitters] == 1:4
            # One localization: MAT reads each column as a scalar.
            smld = load_smite_2d(write_smd(dir, "one.mat", smite_smd(1)))
            @test length(smld.emitters) == 1
            @test smld.emitters[1].x ≈ (10 - 0.5) * 0.1
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
