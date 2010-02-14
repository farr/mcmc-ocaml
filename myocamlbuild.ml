open Ocamlbuild_plugin

let oUnit_dir = "/Users/farr/Documents/code/oUnit"

let _ = dispatch begin function 
  | After_rules -> 
      ocaml_lib "denest";
      ocaml_lib ~extern:true ~dir:oUnit_dir "oUnit"
  | _ -> ()
end
