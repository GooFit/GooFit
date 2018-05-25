# You should run this as:
# cmake -DINFILE="input.h" -DOUTFILE="output.h" -P make_doc_headers.cmake

file(READ "${INFILE}" INPUT_STR)
string(REGEX MATCH "/\\*\\*.*\\*\\*/[ \t\n]*class [a-zA-Z_]+" INPUT_STR "${INPUT_STR}")
string(REGEX REPLACE "/\\*\\*.*\\*\\*/[ \t\n]*class ([a-zA-Z_]+)" "\\1" PDFNAME "${INPUT_STR}")
string(REGEX REPLACE "/\\*\\*(.*)\\*\\*/[ \t\n]*class [a-zA-Z_]+" "\\1" INPUT_STR "${INPUT_STR}")
string(REPLACE [=[\f$]=] "$" INPUT_STR "${INPUT_STR}")
string(REPLACE [=[\f[]=] "$$" INPUT_STR "${INPUT_STR}")
string(REPLACE [=[\f]]=] "$$" INPUT_STR "${INPUT_STR}")
string(REGEX REPLACE "\\\\class [^\n]*" "" INPUT_STR "${INPUT_STR}")
string(CONCAT INPUT_STR  "#pragma once\n\n"
                         "#include <string>\n\n"
                         "const std::string ${PDFNAME}_docs = R\"raw("
                         "${INPUT_STR}"
                         ")raw\";\n")
                     
file(WRITE "${OUTFILE}" "${INPUT_STR}")
