using Test,JLD2,TOML
include(joinpath(@__DIR__,"..","sample_ingestion_v2","export_identities_native.jl"))
function metadata(;label="AVAL",cls="AVA",confidence=5,lr="L",dv="undefined")
 Dict{String,Any}("label"=>label,"neuron_class"=>cls,"confidence"=>confidence,"LR"=>lr,"DV"=>dv,"region"=>"","roi_id"=>[1])
end
function views(items;keys_override=nothing)
 byrow=Dict{Any,Any}();byclass=Dict{Any,Any}()
 for (i,m) in enumerate(items)
  byrow[i]=m;k=isnothing(keys_override) ? m["neuron_class"] : keys_override[i]
  push!(get!(byclass,k,Any[]),(m,i,0.0))
 end
 (byrow,byclass)
end
function fixture(path,root)
 jldsave(path;dict_match_dict=root)
end
@testset "Native exporter function scope and scientific guards" begin
 mktempdir() do d
  source=joinpath(d,"good.jld2");fixture(source,Dict("synthetic_A"=>views([metadata()]),"synthetic_B"=>views([metadata(label="RIS",cls="RIS",lr="undefined")])))
  output=joinpath(d,"cli_output");script=joinpath(@__DIR__,"..","sample_ingestion_v2","export_identities_native.jl");log=joinpath(d,"cli.log")
  cmd=`$(Base.julia_cmd()) --startup-file=no --project=$(dirname(Base.active_project())) $script $source $output`
  open(log,"w") do io;run(pipeline(cmd,stdout=io,stderr=io));end
  text=read(log,String)
  @test !occursin("soft scope",text)
  @test !occursin("UndefVarError",text)
  receipt=TOML.parsefile(joinpath(output,"native_receipt_PRIVATE.toml"))
  @test receipt["animal_dictionary_keys"]==2
  @test receipt["entries_cross_checked"]==2
  @test isfile(joinpath(output,"labels_PRIVATE.tsv"))
  @test_throws ErrorException export_identities(source,output)
  @test lookup_relation("ASEL",metadata(label="ASEL",cls="ASE"))=="individual_label"
  @test lookup_relation("CEPD",metadata(label="CEPDL",cls="CEP",dv="D"))=="dorsoventral_refinement"
  @test_throws ErrorException lookup_relation("UNEXPLAINED",metadata())
  bad=joinpath(d,"bad_conf.jld2");fixture(bad,Dict("synthetic"=>views([metadata(confidence=6)])))
  r=export_identities(bad,joinpath(d,"bad_conf_output"))
  @test r["quarantined_rows"]==1
  @test r["valid_identity_rows"]==0
  @test r["quarantine_reasons"]["invalid_confidence"]==1
  @test occursin("\t6\t",read(joinpath(d,"bad_conf_output","labels_PRIVATE.tsv"),String))
  duplicate=joinpath(d,"duplicate.jld2");fixture(duplicate,Dict("synthetic"=>views([metadata(),metadata(confidence=3)])))
  r=export_identities(duplicate,joinpath(d,"duplicate_output"))
  @test r["duplicate_individual_label_groups"]==1
  @test r["quarantined_rows"]==2
  @test r["valid_identity_rows"]==0
  inconsistent=views([metadata()]);inconsistent[2]["AVA"][1]=(metadata(confidence=3),1,0.)
  conflict=joinpath(d,"conflict.jld2");fixture(conflict,Dict("synthetic"=>inconsistent))
  @test_throws ErrorException export_identities(conflict,joinpath(d,"conflict_output"))
  @test !ispath(joinpath(d,"conflict_output"))
  repeated=views([metadata()]);push!(repeated[2]["AVA"],repeated[2]["AVA"][1])
  repeatfile=joinpath(d,"repeated.jld2");fixture(repeatfile,Dict("synthetic"=>repeated))
  @test_throws ErrorException export_identities(repeatfile,joinpath(d,"repeat_output"))
  lr=joinpath(d,"lr.jld2");fixture(lr,Dict("synthetic"=>views([metadata(lr="R")])))
  r=export_identities(lr,joinpath(d,"lr_output"))
  @test r["quarantine_reasons"]["inconsistent_LR"]==1
 end
end
