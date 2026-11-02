using JLD2,TOML
function inspect_schema(input,output)
    root=jldopen(input,"r") do f; f["dict_match_dict"]; end
    mismatches=Dict{String,Any}[];missing=Dict{String,Any}[];scales=Dict{String,Int}();fieldnames=Set{String}();rows=0;entries=0;conflicts=0;ambiguous_rows=0
    animals=String[]
    for (animal,views) in root
        push!(animals,string(animal));byrow,byclass=views
        counts=Dict{Int,Int}()
        for (key,list) in byclass
            for entry in list
                meta,idx=entry[1],entry[2];entries+=1;counts[idx]=get(counts,idx,0)+1
                if !haskey(byrow,idx);push!(missing,Dict("animal"=>string(animal),"row"=>idx));continue;end
                for k in ["label","neuron_class","confidence","LR","DV","region"]
                    if !haskey(meta,k)||!haskey(byrow[idx],k)||!isequal(meta[k],byrow[idx][k]);conflicts+=1;end
                end
                if string(meta["neuron_class"])!=string(key)
                    push!(mismatches,Dict("animal"=>string(animal),"row"=>idx,"class_lookup_key"=>string(key),"metadata_neuron_class"=>string(meta["neuron_class"]),"label"=>string(meta["label"]),"LR"=>string(meta["LR"]),"confidence"=>string(meta["confidence"])))
                end
            end
        end
        ambiguous_rows+=count(x->x>1,values(counts))
        for (idx,meta) in byrow
            rows+=1;union!(fieldnames,string.(collect(keys(meta))));s=string(meta["confidence"]);scales[s]=get(scales,s,0)+1
        end
    end
    result=Dict("dictionary_animal_keys"=>length(root),"animal_keys_PRIVATE"=>sort(animals),"row_entries"=>rows,"class_view_entries"=>entries,"class_key_mismatch_count"=>length(mismatches),"mismatches_PRIVATE"=>mismatches,"missing_rows_PRIVATE"=>missing,"two_view_field_conflicts"=>conflicts,"repeated_class_view_rows"=>ambiguous_rows,"confidence_counts"=>scales,"metadata_field_names"=>sort(collect(fieldnames)))
    open(output,"w") do io;TOML.print(io,result);end
    println("Native schema audit: animal keys=",length(root),"; row entries=",rows,"; class-key mismatches=",length(mismatches),"; view conflicts=",conflicts)
end
inspect_schema(ARGS[1],ARGS[2])
