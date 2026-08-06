#!/bin/bash

# Usage: ./compile_simple.sh monfichier.tex
# Simple compilation for documents with bibliography

if [ $# -lt 1 ]; then
  echo "Usage: $0 fichier.tex"
  exit 1
fi

TEXFILE="$1"
BASENAME="${TEXFILE%.tex}"

# Compilation LaTeX + bibliography -> PDF (using pdflatex + bibtex)
pdflatex "$TEXFILE"
bibtex "$BASENAME"
pdflatex "$TEXFILE"
pdflatex "$TEXFILE"  # Final pass for references and labels

echo "Compilation terminée : ${BASENAME}.pdf"

# Clean up LaTeX temporary files generated during compilation
rm -f \
  "${BASENAME}.aux" \
  "${BASENAME}.bbl" \
  "${BASENAME}.blg" \
  "${BASENAME}.lof" \
  "${BASENAME}.lot" \
  "${BASENAME}.log" \
  "${BASENAME}.out" \
  "${BASENAME}.spl" \
  "${BASENAME}.toc" \
  "${BASENAME}.synctex.gz"
