using JLD2,TOML
function main(input,output)
 root=jldopen(input,"r") do f;f["dict_match_dict"];end
 unknown=[];duplicates=[];lr_conflicts=[];invalid=[];ambiguous=[];counts=Dict("coarse"=>0,"individual"=>0,"dv_refinement"=>0,"other"=>0)
 for (animal,views) in root
  byrow,byclass=views;seen=Dict{String,Vector{Int}}()
  for (key,entries) in byclass, entry in entries
   m,idx=entry[1],entry[2];k=string(key);c=string(m["neuron_class"]);label=string(m["label"]);dv=string(m["DV"]);lr=string(m["LR"])
   method=k==c ? "coarse" : k==label ? "individual" : dv in ["D","V"] && k==c*dv ? "dv_refinement" : "other"
   counts[method]+=1
   method=="other" && push!(unknown,Dict("animal"=>string(animal),"row"=>idx,"key"=>k,"class"=>c,"label"=>label,"DV"=>dv,"LR"=>lr))
   if lr in ["L","R"] && !endswith(label,lr);push!(lr_conflicts,Dict("animal"=>string(animal),"row"=>idx,"label"=>label,"LR"=>lr));end
   if !(m["confidence"] isa Real && isfinite(m["confidence"]) && 1<=m["confidence"]<=5);push!(invalid,Dict("animal"=>string(animal),"row"=>idx,"confidence"=>string(m["confidence"])));end
   push!(get!(seen,label,Int[]),idx)
  end
  for (label,rs) in seen
   if length(rs)>1
    push!(duplicates,Dict("animal"=>string(animal),"label"=>label,"rows"=>rs,"confidence"=>[string(byrow[r]["confidence"]) for r in rs]))
   end
  end
 end
 d=Dict("class_refinement_categories"=>counts,"unknown_key_relations_PRIVATE"=>unknown,"duplicate_labels_PRIVATE"=>duplicates,"LR_conflicts_PRIVATE"=>lr_conflicts,"invalid_confidence_PRIVATE"=>invalid)
 open(output,"w") do io;TOML.print(io,d);end
 println("Key relations=",counts,"; unknown=",length(unknown),"; duplicate labels=",length(duplicates),"; LR conflicts=",length(lr_conflicts),"; invalid confidence=",length(invalid))
end
main(ARGS[1],ARGS[2])
