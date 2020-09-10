# You should run this as:
# cmake -DINFILE="input.h" -DOUTFILE="output.h" -P make_doc_headers.cmake

set(OUTPUT_STR "#pragma once\n\n#include <string>\n\n")

file(READ "${INFILE}" INPUT_FULL_STR)

string(REPLACE ";" "_CMAKE_SIMICOLON_TOKEN_" INPUT_FULL_STR "${INPUT_FULL_STR}")

# Class docs
string(REGEX MATCHALL "/\\*\\*.*\\*\\*/[ \t\n]*(class|struct) [a-zA-Z0-9_]+" INPUT_LIST
             "${INPUT_FULL_STR}")

foreach(INPUT_STR IN LISTS INPUT_LIST)
  string(REPLACE "_CMAKE_SIMICOLON_TOKEN_" ";" INPUT_STR "${INPUT_STR}")
  string(REGEX REPLACE "/\\*\\*.*\\*\\*/[ \t\n]*(class|struct) ([a-zA-Z0-9_]+)" "\\2" PDFNAME
                       "${INPUT_STR}")
  string(REGEX REPLACE "/\\*\\*(.*)\\*\\*/[ \t\n]*(class|struct) [a-zA-Z0-9_]+" "\\1" INPUT_STR
                       "${INPUT_STR}")

  # Inline math mode
  string(REPLACE [=[\f$]=] "$" INPUT_STR "${INPUT_STR}")

  # Centered math mode
  string(REPLACE [=[\f[]=] "$$" INPUT_STR "${INPUT_STR}")
  string(REPLACE [=[\f]]=] "$$" INPUT_STR "${INPUT_STR}")

  # Footnotes
  # \anchor footnote1 1 (: is in both)
  string(REGEX REPLACE [=[\\anchor footnote[0-9]+ ([0-9]+)]=] "[^\\1]" INPUT_STR "${INPUT_STR}")
  # (\ref footnote2 "2")
  string(REGEX REPLACE [=[ ?\(\\ref footnote[0-9]+ "([0-9]+)"\)]=] "[^\\1]" INPUT_STR
                       "${INPUT_STR}")

  string(CONCAT OUTPUT_STR "${OUTPUT_STR}" "const std::string ${PDFNAME}_docs = R\"raw("
                "${INPUT_STR}" ")raw\";\n")

endforeach()

file(WRITE "${OUTFILE}" "${OUTPUT_STR}")
