# Native author loading semantics; no Python traversal of serialized Julia objects.
using JLD2
using SHA
using TOML

function safe(value)
    text=string(value)
    occursin(r"[\t\r\n]",text) && error("Unexpected metadata delimiter")
    return text
end
function field(metadata,key)
    haskey(metadata,key) || error("Missing published metadata field: $key")
    return metadata[key]
end
function lookup_relation(key,metadata)
    k=string(key);c=string(field(metadata,"neuron_class"));label=string(field(metadata,"label"));dv=string(field(metadata,"DV"))
    k==c && return "coarse_class"
    k==label && return "individual_label"
    dv in ("D","V") && k==c*dv && return "dorsoventral_refinement"
    error("Unknown class-lookup relationship; do not reinterpret metadata")
end

function export_identities(source::AbstractString,destination::AbstractString)
    source=abspath(source);destination=abspath(destination)
    isfile(source) || error("Original identity dictionary absent")
    ispath(destination) && error("Output exists; preserve it and use a versioned directory")
    labels=jldopen(source,"r") do file
        file["dict_match_dict"]
    end
    labels isa AbstractDict || error("Author root is not a dictionary")
    lines=String["animal\ttrace_row_1based\tlabel\tneuron_class\tconfidence\tLR\tDV\tregion\tclass_lookup_key\tlookup_relation\tidentity_valid\trejection_reason"]
    checked=0;animals=0;quarantined=0;duplicate_groups=0;duplicate_rows=0;valid_rows=0
    metadata_fields=Set{String}();confidence_counts=Dict{String,Int}();relations=Dict{String,Int}();reasons_count=Dict{String,Int}()
    for animal in sort(collect(keys(labels)),by=string)
        views=labels[animal]
        length(views)==2 || error("Expected the two author lookup views")
        byrow,byclass=views[1],views[2]
        byrow isa AbstractDict && byclass isa AbstractDict || error("Unexpected native dictionary views")
        row_keys=Dict{Int,Tuple{String,String}}()
        for (key,entries) in byclass
            for entry in entries
                length(entry)==3 || error("Unexpected class lookup tuple")
                metadata,index=entry[1],entry[2]
                index isa Integer && index>0 || error("Trace-row index is not a positive integer")
                haskey(byrow,index) || error("Class view row absent from row view")
                haskey(row_keys,Int(index)) && error("Repeated or ambiguous trace-row class lookup")
                for name in ("label","neuron_class","confidence","LR","DV","region")
                    isequal(field(metadata,name),field(byrow[index],name)) || error("Conflicting row/class metadata: $name")
                end
                relation=lookup_relation(key,metadata)
                row_keys[Int(index)]=(string(key),relation)
                relations[relation]=get(relations,relation,0)+1
                checked+=1
            end
        end
        label_rows=Dict{String,Vector{Int}}()
        for (index,metadata) in byrow
            index isa Integer && index>0 || error("Invalid native trace index")
            haskey(row_keys,Int(index)) || error("Row entry absent from class view")
            push!(get!(label_rows,string(field(metadata,"label")),Int[]),Int(index))
        end
        repeated=Set{Int}()
        for indices in values(label_rows)
            if length(indices)>1
                duplicate_groups+=1;duplicate_rows+=length(indices);union!(repeated,indices)
            end
        end
        for index in sort(collect(keys(byrow)))
            metadata=byrow[index];union!(metadata_fields,string.(collect(keys(metadata))))
            label=string(field(metadata,"label"));cls=string(field(metadata,"neuron_class"));lr=string(field(metadata,"LR"));dv=string(field(metadata,"DV"));confidence=field(metadata,"confidence")
            score=string(confidence);confidence_counts[score]=get(confidence_counts,score,0)+1
            reasons=String[]
            confidence isa Real && isfinite(confidence) && isinteger(confidence) && 1<=confidence<=5 || push!(reasons,"invalid_confidence")
            index in repeated && push!(reasons,"duplicate_individual_label")
            isempty(label) && push!(reasons,"missing_individual_label")
            isempty(cls) && push!(reasons,"missing_class")
            lr in ("L","R","undefined") || push!(reasons,"unknown_LR")
            dv in ("D","V","undefined") || push!(reasons,"unknown_DV")
            lr in ("L","R") && !endswith(label,lr) && push!(reasons,"inconsistent_LR")
            if dv in ("D","V")
                expected=cls*dv*(lr in ("L","R") ? lr : "")
                label==expected || push!(reasons,"inconsistent_DV")
            end
            for reason in reasons;reasons_count[reason]=get(reasons_count,reason,0)+1;end
            valid=isempty(reasons)
            if valid;valid_rows+=1;else;quarantined+=1;end
            key,relation=row_keys[Int(index)]
            push!(lines,join(safe.([animal,index,label,cls,confidence,lr,dv,field(metadata,"region"),key,relation,valid,join(reasons,';')]),'\t'))
        end
        animals+=1
    end
    # Structural conflicts abort; invalid/ambiguous records are losslessly exported
    # with an explicit exclusion flag. They are never clamped or counted as valid.
    mkpath(destination)
    table=joinpath(destination,"labels_PRIVATE.tsv")
    open(table,"w") do io;write(io,join(lines,"\n")*"\n");end
    receipt=Dict("status"=>(quarantined==0 ? "NATIVE_VALIDATION_PASS" : "NATIVE_STRUCTURE_PASS_WITH_QUARANTINED_RECORDS"),"source_sha256"=>bytes2hex(open(sha256,source)),"julia_version"=>string(VERSION),"JLD2_version"=>string(Base.pkgversion(JLD2)),"animal_dictionary_keys"=>animals,"entries_cross_checked"=>checked,"rows_exported"=>length(lines)-1,"valid_identity_rows"=>valid_rows,"quarantined_rows"=>quarantined,"duplicate_individual_label_groups"=>duplicate_groups,"duplicate_individual_label_rows"=>duplicate_rows,"quarantine_reasons"=>reasons_count,"confidence_counts"=>confidence_counts,"class_lookup_relations"=>relations,"metadata_field_names"=>sort(collect(metadata_fields)),"labels_sha256"=>bytes2hex(open(sha256,table)))
    open(joinpath(destination,"native_receipt_PRIVATE.toml"),"w") do io;TOML.print(io,receipt);end
    println("Native identity views agree: ",animals," animals; ",checked," rows checked; ",quarantined," rows explicitly quarantined.")
    return receipt
end
function main(args=ARGS)
    length(args)==2 || error("Provide original JLD2 and a new private output directory")
    return export_identities(args[1],args[2])
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main()
end
