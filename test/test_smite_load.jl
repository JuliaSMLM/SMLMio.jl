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
            # A single-precision PixelSize matches the decimal it was written from.
            smd = write_smd(dir, "px32.mat", smite_smd(3; pixel_size=0.0978f0))
            @test load_smite_2d(smd; pixel_size=0.0978).camera.pixel_edges_x[2] ≈ 0.0978 rtol = 1e-6
            @test_throws r"but pixel_size = 0.1" load_smite_2d(smd; pixel_size=0.1)
        end
    end

    @testset "load_smite_2d on a 3D file says it leaves out Z" begin
        mktempdir() do dir
            smd = write_smd(dir, "three.mat", smite_smd_3d(3))
            @test_logs (:warn, r"NDims 3.*load_smite_3d") load_smite_2d(smd)
            smld = quiet(() -> load_smite_2d(smd))
            @test smld.emitters[1] isa Emitter2DFit
            @test haskey(smld.metadata, "Z")
            @test_logs load_smite_2d(write_smd(dir, "two.mat", smite_smd(3)))   # 2D: no warning
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

    @testset "complex values" begin
        mktempdir() do dir
            # A localization with a non-zero imaginary part in a position, photon or background
            # field is dropped with a warning; a zero imaginary part loads as its real part.
            s = smite_smd(5)
            s["X"] = col(ComplexF32[10, 17.5, 25 + 0im, 32.5, 40])
            s["Photons_SE"] = col(ComplexF32[30, 30, 30 + 2im, 30, 30])
            s["Bg"] = col(ComplexF32[5, 5, 5, 5, 5 - 1im])
            smd = write_smd(dir, "complex.mat", s)
            @test_logs (:warn, r"Found 2 emitters with non-zero imaginary") load_smite_2d(smd)
            smld = quiet(() -> load_smite_2d(smd))
            @test [e.id for e in smld.emitters] == [1, 2, 4]
            @test [e.x for e in smld.emitters] ≈ (Float64[10, 17.5, 32.5] .- 0.5) .* 0.1
            @test smld.metadata["complex_fields_removed"] == true
            @test sort(smld.metadata["complex_fields"]) == ["Bg", "Photons_SE"]
            @test (smld.metadata["original_emitter_count"], smld.metadata["removed_emitter_count"]) == (5, 2)
            # Complex values in a field loaded as its default (one per dataset) drop nothing.
            s = smite_smd(3); s["NDatasets"] = 2.0
            s["Bg_SE"] = col(ComplexF32[0.5 + 1im, 0.5])
            smld = quiet(() -> load_smite_2d(write_smd(dir, "complexdefault.mat", s)))
            @test length(smld.emitters) == 3
            @test !haskey(smld.metadata, "complex_fields_removed")
            # Row vectors: the dropped localization is the right one.
            s = smite_smd(4)
            s["Y"] = col(ComplexF32[5, 20, 35 + 1im, 50])
            for (k, v) in s
                v isa AbstractMatrix && size(v, 2) == 1 && size(v, 1) == 4 && (s[k] = permutedims(v))
            end
            smld = quiet(() -> load_smite_2d(write_smd(dir, "complexrows.mat", s)))
            @test [e.id for e in smld.emitters] == [1, 2, 4]
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

    @testset "save_smite then load: the same data" begin
        # SMITE has no emitter id; the loaders number emitters by row, so `skip` drops `id` when the
        # data were not loaded from SMITE.
        same_emitters(a, b; skip=()) = length(a) == length(b) && all(zip(a, b)) do (p, q)
            all(f -> f in skip || isequal(getfield(p, f), getfield(q, f)) ||
                     (getfield(p, f) isa AbstractFloat && getfield(p, f) ≈ getfield(q, f)),
                fieldnames(typeof(p)))
        end
        mktempdir() do dir
            # load, save, load (2D, 3D without ZSize, a file with defaulted fields)
            s2 = smite_smd(4; pixel_size=0.108)
            s3 = smite_smd_3d(4; pixel_size=0.0978)
            sd = smite_smd(4); sd["ConnectID"] = zeros(UInt32, 0, 0); sd["Bg_SE"] = zeros(0, 0)
            for (name, s, load) in (("2d", s2, load_smite_2d), ("3d", s3, load_smite_3d), ("def", sd, load_smite_2d))
                a = load(write_smd(dir, "$name.mat", s))
                save_smite(a, dir, "$name-saved.mat")
                b = load(SmiteSMD(dir, "$name-saved.mat"))
                @test same_emitters(a.emitters, b.emitters)
                @test b.camera.pixel_edges_x ≈ a.camera.pixel_edges_x
                @test b.camera.pixel_edges_y ≈ a.camera.pixel_edges_y
                @test (b.n_frames, b.n_datasets) == (a.n_frames, a.n_datasets)
                saved = matread(joinpath(dir, "$name-saved.mat"))["SMD"]
                @test (saved["XSize"], saved["YSize"], saved["PixelSize"]) == (64, 64, s["PixelSize"])
                @test saved["NDims"] == (name == "3d" ? 3 : 2)
                @test !haskey(saved, "data_size") && !haskey(saved, "pixel_size")
                @test !haskey(saved, "smite_defaulted_fields")
                @test vec(saved["DriftX"]) == vec(s["DriftX"])   # other SMITE fields kept
            end
            # Any SMLD saves, not only a SmiteSMLD, and loads back.
            cam = IdealCamera(1:32, 1:16, 0.1)
            es = [Emitter2DFit{Float64}(1.0, 0.5, 1000.0, 5.0, 0.02, 0.03, 30.0, 0.5; frame=2, dataset=1),
                  Emitter2DFit{Float64}(2.5, 1.2, 900.0, 4.0, 0.02, 0.03, 30.0, 0.5; frame=7, dataset=1)]
            basic = BasicSMLD(es, cam, 10, 1, Dict{String,Any}())
            save_smite(basic, dir, "basic.mat")
            b = load_smite_2d(SmiteSMD(dir, "basic.mat"))
            @test same_emitters(es, b.emitters; skip=(:id,))
            @test [e.id for e in b.emitters] == [1, 2]
            @test b.camera.pixel_edges_x ≈ cam.pixel_edges_x
            @test b.camera.pixel_edges_y ≈ cam.pixel_edges_y
            # SMITE has one square pixel size: a camera with other pixels is an error.
            rect = BasicSMLD(es, IdealCamera(32, 16, (0.1, 0.2)), 10, 1, Dict{String,Any}())
            @test_throws r"one square pixel size" save_smite(rect, dir, "rect.mat")
        end
    end

    @testset "save_smite keeps SMITE's other fields with their localizations" begin
        mktempdir() do dir
            # Per-localization fields kept in metadata (PValue, LogLikelihood, ...) follow the
            # emitters through a filter and through complex screening, by the file row in `id`.
            s = smite_smd(5)
            s["PValue"] = col(Float32[0.1, 0.2, 0.3, 0.4, 0.5])
            s["LogLikelihood"] = col([-1.0, -2.0, -3.0, -4.0, -5.0])
            a = load_smite_2d(write_smd(dir, "a.mat", s))
            @test sort(a.metadata["smite_per_localization_fields"]) == ["LogLikelihood", "PValue", "ThreshFlag"]
            kept = SMLMCore.@filter(a, frame >= 3)        # rows 3, 4, 5
            save_smite(kept, dir, "kept.mat")
            k = matread(joinpath(dir, "kept.mat"))["SMD"]
            @test vec(k["PValue"]) == Float32[0.3, 0.4, 0.5]
            @test vec(k["LogLikelihood"]) == [-3.0, -4.0, -5.0]
            @test length(vec(k["ThreshFlag"])) == 3
            @test vec(k["X"]) ≈ vec(s["X"])[3:5]
            # Complex screening drops row 2; the saved fields drop it too.
            c = copy(s)
            c["X"] = col(ComplexF32[10, 17.5 + 1im, 25, 32.5, 40])
            b = quiet(() -> load_smite_2d(write_smd(dir, "complex.mat", c)))
            @test [e.id for e in b.emitters] == [1, 3, 4, 5]
            save_smite(b, dir, "screened.mat")
            sc = matread(joinpath(dir, "screened.mat"))["SMD"]
            @test vec(sc["PValue"]) == Float32[0.1, 0.3, 0.4, 0.5]
            # Data sets combined from several files cannot be traced to one file's rows: the fields
            # are dropped with a warning, never written next to the wrong localizations.
            both = SMLMCore.cat_smld(a, kept)
            @test_logs (:warn, r"not saved") match_mode = :any save_smite(both, dir, "both.mat")
            bo = matread(joinpath(dir, "both.mat"))["SMD"]
            @test !haskey(bo, "PValue") && !haskey(bo, "LogLikelihood")
            @test length(vec(bo["X"])) == 8
        end
    end

    @testset "save_smite skips the loaders' own metadata" begin
        mktempdir() do dir
            # Loading B must describe B, not the file B was saved from.
            c = smite_smd(3)
            c["X"] = col(ComplexF32[10, 25 + 1im, 40])
            a = quiet(() -> load_smite_2d(write_smd(dir, "A.mat", c)))
            @test a.metadata["complex_fields_removed"] == true
            save_smite(a, dir, "B.mat")
            saved = matread(joinpath(dir, "B.mat"))["SMD"]
            for k in ("original_file", "complex_fields_removed", "complex_fields", "original_emitter_count",
                      "removed_emitter_count", "smite_per_localization_fields")
                @test !haskey(saved, k)
            end
            b = load_smite_2d(SmiteSMD(dir, "B.mat"))
            @test b.metadata["original_file"] == "B.mat"
            @test !haskey(b.metadata, "complex_fields_removed")
            @test !occursin("Complex", sprint(show, MIME("text/plain"), b))
            # A file field that shares a loader key's name does not replace the loader's value.
            s = smite_smd(2); s["original_file"] = "elsewhere.mat"
            @test load_smite_2d(write_smd(dir, "C.mat", s)).metadata["original_file"] == "C.mat"
        end
    end

    @testset "save_smite writes untouched defaulted fields back as SMITE had them" begin
        mktempdir() do dir
            # SMITE treats an empty ConnectID (no frame connection) differently from zeros.
            s = smite_smd(3); s["ConnectID"] = zeros(UInt32, 0, 0)
            s["NDatasets"] = 2.0; s["Photons_SE"] = col(Float32[31, 32])
            a = quiet(() -> load_smite_2d(write_smd(dir, "a.mat", s)))
            save_smite(a, dir, "b.mat")
            saved = matread(joinpath(dir, "b.mat"))["SMD"]
            @test isempty(saved["ConnectID"])
            @test vec(saved["Photons_SE"]) == Float32[31, 32]
            b = quiet(() -> load_smite_2d(SmiteSMD(dir, "b.mat")))
            @test b.metadata["smite_defaulted_fields"] == a.metadata["smite_defaulted_fields"]
            # Values set since loading are written.
            a.emitters[1].track_id = 7
            save_smite(a, dir, "c.mat")
            @test vec(matread(joinpath(dir, "c.mat"))["SMD"]["ConnectID"]) == [7.0, 0.0, 0.0]
        end
    end

    @testset "save_smite cameras and number types" begin
        mktempdir() do dir
            # Float32 cameras, whose edges carry Float32 rounding, save and load back.
            e32 = [Emitter2DFit{Float32}(1.0f0, 2.0f0, 1000.0f0, 5.0f0, 0.02f0, 0.03f0, 30.0f0, 0.5f0),
                   Emitter2DFit{Float32}(40.0f0, 25.5f0, 900.0f0, 4.0f0, 0.02f0, 0.03f0, 30.0f0, 0.5f0)]
            for cam in (IdealCamera(1:512, 1:512, 0.1f0), SMLMCore.SCMOSCamera(256, 256, 0.1f0, 1.6f0))
                save_smite(BasicSMLD(e32, cam, 1, 1, Dict{String,Any}()), dir, "f32.mat")
                b = load_smite_2d(SmiteSMD(dir, "f32.mat"))
                @test [e.x for e in b.emitters] ≈ [1.0, 40.0] rtol = 1e-6
                @test [e.y for e in b.emitters] ≈ [2.0, 25.5] rtol = 1e-6
                @test length(b.camera.pixel_edges_x) == length(cam.pixel_edges_x)
            end
            # Float32 emitters keep Float32 position columns (MATLAB single), as before.
            saved = matread(joinpath(dir, "f32.mat"))["SMD"]
            @test eltype(saved["X"]) == Float32 && eltype(saved["X_SE"]) == Float32
            # The saved PixelSize is the camera's pixel, not a rounded mean over its width.
            e64 = [Emitter2DFit{Float64}(1.0, 2.0, 1000.0, 5.0, 0.02, 0.03, 30.0, 0.5)]
            save_smite(BasicSMLD(e64, IdealCamera(1:200, 1:200, 0.0954), 1, 1, Dict{String,Any}()), dir, "px.mat")
            saved = matread(joinpath(dir, "px.mat"))["SMD"]
            @test saved["PixelSize"] === 0.0954
            # Counts and indices are doubles, as SMITE writes them.
            @test all(k -> saved[k] isa Float64, ("XSize", "YSize", "NDims", "NFrames", "NDatasets"))
            @test all(k -> eltype(saved[k]) == Float64, ("FrameNum", "DatasetNum", "ConnectID"))
            # SMITE has no origin: a camera whose first edge is not at 0 warns, and the file loads
            # back with the camera at 0 and the positions moved with it.
            off = IdealCamera(collect(1.0:0.1:7.4), collect(2.0:0.1:8.4))
            e = [Emitter2DFit{Float64}(1.05, 2.25, 1000.0, 5.0, 0.02, 0.03, 30.0, 0.5)]
            @test_logs (:warn, r"first pixel edge is at \(1.0, 2.0\)") save_smite(BasicSMLD(e, off, 1, 1, Dict{String,Any}()), dir, "off.mat")
            b = load_smite_2d(SmiteSMD(dir, "off.mat"))
            @test b.emitters[1].x ≈ 0.05 && b.emitters[1].y ≈ 0.25
            @test b.camera.pixel_edges_x[1] == 0
        end
    end
end
