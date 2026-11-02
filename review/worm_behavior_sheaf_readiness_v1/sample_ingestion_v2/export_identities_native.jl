# Author loading semantics: f=jldopen(dict_path); NeuroPAL_labels=f["dict_match_dict"]; close(f).
# No Python traversal of Julia objects. No physiological modelling.
using JLD2
using SHA
using TOML
length(ARGS)==2 || error("Provide original JLD2 and a new private output directory")
source=abspath(ARGS[1]); destination=abspath(ARGS[2])
isfile(source) || error("Original identity dictionary absent")
ispath(destination) && error("Output exists; preserve it and use a new versioned directory")
f=jldopen(source,"r")
NeuroPAL_labels=try
    f["dict_match_dict"]
finally
    close(f)
end
NeuroPAL_labels isa AbstractDict || error("Author root is not a dictionary")
function safe(v)
    s=string(v)
    occursin(r"[\t\r\n]",s) && error("Unexpected delimiter; no silent escaping")
    s
end
function field(meta,key)
    haskey(meta,key) || error("Missing published metadata field: $key")
    meta[key]
end
lines=String["animal\ttrace_row_1based\tlabel\tneuron_class\tconfidence\tLR\tDV\tregion"]
checked=0; animals=0; duplicates=0
metadata_fields=Set{String}()
for animal in sort(collect(keys(NeuroPAL_labels)),by=string)
    views=NeuroPAL_labels[animal]
    length(views)==2 || error("Expected the two author dictionary views")
    byrow,byclass=views[1],views[2]
    byrow isa AbstractDict && byclass isa AbstractDict || error("Unexpected native dictionary views")
    cross=Dict{Int,Vector{Any}}()
    for (cls,entries) in byclass
        for entry in entries
            length(entry)==3 || error("Unexpected class lookup tuple")
            meta,idx=entry[1],entry[2]
            idx isa Integer || error("Trace-row key is not an integer")
            haskey(byrow,idx) || error("Class view row absent from row view")
            for k in ["label","neuron_class","confidence","LR","DV","region"]
                isequal(field(meta,k),field(byrow[idx],k)) || error("Conflicting two-view metadata at $animal row $idx: $k")
            end
            string(field(meta,"neuron_class"))==string(cls) || error("Class-key mismatch")
            push!(get!(cross,Int(idx),Any[]),meta)
            checked+=1
        end
    end
    for idx in sort(collect(keys(byrow)))
        idx isa Integer && idx>0 || error("Invalid native trace index")
        haskey(cross,Int(idx)) || error("Row view entry absent from class view")
        length(cross[Int(idx)])==1 || error("Repeated/ambiguous class lookup for a trace row")
        meta=byrow[idx]
        union!(metadata_fields,string.(collect(keys(meta))))
        confidence=field(meta,"confidence")
        confidence isa Real && isfinite(confidence) && 1<=confidence<=5 || error("Unexpected identity-confidence scale")
        push!(lines,join(safe.([animal,idx,field(meta,"label"),field(meta,"neuron_class"),confidence,field(meta,"LR"),field(meta,"DV"),field(meta,"region")]),'\t'))
    end
    animals+=1
end
# Write only after complete native validation; individual metadata stays private.
mkpath(destination)
open(joinpath(destination,"labels_PRIVATE.tsv"),"w") do io
    write(io,join(lines,"\n")*"\n")
end
receipt=Dict("status"=>"NATIVE_TWO_VIEW_VALIDATION_PASS","source_sha256"=>bytes2hex(open(sha256,source)),"julia_version"=>string(VERSION),"JLD2_version"=>string(Base.pkgversion(JLD2)),"animal_dictionary_keys"=>animals,"entries_cross_checked"=>checked,"rows_exported"=>length(lines)-1,"metadata_field_names"=>sort(collect(metadata_fields)),"labels_sha256"=>bytes2hex(open(sha256,joinpath(destination,"labels_PRIVATE.tsv"))))
open(joinpath(destination,"native_receipt_PRIVATE.toml"),"w") do io
    TOML.print(io,receipt)
end
println("Native Julia identity views agree; exported private metadata for ",animals," dictionary keys. No model fitted.")
