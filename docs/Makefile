# Minimal Makefile for Sphinx documentation
SPHINXOPTS   ?=
SPHINXBUILD  ?= sphinx-build
SOURCEDIR    = source
BUILDDIR     = build

# @NOTE(ds): Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help latexpdf Makefile

# latexpdf: run Sphinx latexpdf then a second pdflatex pass so that the
# List of Figures is populated (lof entries are written on pass 1, read on pass 2).
latexpdf:
	@$(SPHINXBUILD) -M latexpdf "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	@cd "$(BUILDDIR)/latex" && pdflatex -interaction=nonstopmode pyfraglib_documentation.tex > /dev/null
	@echo "latexpdf: second pdflatex pass done — List of Figures populated."

# Catch-all target: route all unknown targets to Sphinx-Makefile
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
