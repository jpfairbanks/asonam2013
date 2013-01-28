LATEX = pdflatex
VIEWER = evince
SHELL = /bin/bash
LTXFLAGS =
DOCNAME = writeup
BIBENGINE = bibtex
tex: $(DOCNAME).md
	pandoc --read=markdown --write=latex --output=$(DOCNAME).tex --variable=documentclass:scrartcl --include-in-header=header.tex --include-before-body=includebody.tex --self-contained --table-of-contents $(DOCNAME).md
pdf: $(DOCNAME).tex
	$(LATEX) $(DOCNAME).tex
bib: $(DOCNAME).tex
	$(BIBENGINE) $(DOCNAME)
	$(LATEX) $(DOCNAME).tex
	$(LATEX) $(DOCNAME).tex
view: pdf
	$(VIEWER) $(DOCNAME).pdf
clean:
	rm $(DOCNAME).pdf $(DOCNAME).aux $(DOCNAME).toc $(DOCNAME).html $(DOCNAME).dvi $(DOCNAME).log
