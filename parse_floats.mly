%token <float> FLOAT
%token EOL EOF
%start parse_floats
%type <(float array) list> parse_floats

%%

parse_floats:
  lines ennd { List.rev $1 }
;

ennd:
  EOF { () }
| EOL ennd { () }

lines: 
   { [] }
| lines line { $2 :: $1 }

line:
  floats EOL { Array.of_list (List.rev $1) }

floats:
   { [] }
| floats FLOAT { $2 :: $1 }
