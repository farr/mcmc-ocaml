open Ocamlbuild_plugin

let oUnit_dir = "/Users/farr/Documents/code/oUnit"
let lacaml_dir = "/Users/farr/Documents/code/lacaml/lib"

let _ = dispatch begin function 
  | After_rules -> 
      ocaml_lib "mcmc";
      ocaml_lib ~extern:true ~dir:oUnit_dir "oUnit";
      ocaml_lib ~extern:true ~dir:lacaml_dir "lacaml"        
  | _ -> ()
end
