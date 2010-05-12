{ open Parse_floats }

let digit = ['0'-'9']
let uint = digit+
let pint = '+' uint
let nint = '-' uint
let int = uint | pint | nint

let exponent = 'e' | 'E'
let dot = '.' | ','

let eol = '\n'

let wspace = [' ' '\t' '\r']+

let float = int | (int dot int) | (int dot int exponent int) | (int exponent int)

let inf = "inf"

let neg_inf = "-inf"

let nan = "nan"

rule lex_floats = parse 
    eol { EOL }
  | wspace { lex_floats lexbuf }
  | float as x { FLOAT(float_of_string x) }
  | inf { FLOAT(infinity) }
  | neg_inf { FLOAT(neg_infinity) }
  | nan { FLOAT( 0.0 /. 0.0 ) }
  | eof { EOF }
  | _ { lex_floats lexbuf }
