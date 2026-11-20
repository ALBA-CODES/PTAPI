# PTAPI v2 - relatorio tecnico (../PTAPI-v2.pdf)
#   make         gera o PDF (LuaLaTeX, tres passagens por causa do sumario)
#   make clean   apaga os arquivos intermediarios
#
# Precisa de TeX Live com LuaLaTeX e dos pacotes: babel-portuges, listings,
# booktabs, pgf, tcolorbox, titlesec, enumitem, microtype, unicode-math,
# libertinus-fonts e dejavu. Se o TeX Live estiver instalado na pasta pessoal
# (~/texlive), ele e' encontrado sem alterar o PATH.
TEXBIN := $(firstword $(wildcard $(HOME)/texlive/bin/*))
ifneq ($(TEXBIN),)
export PATH := $(TEXBIN):$(PATH)
endif
B   = build
TEX = lualatex -interaction=nonstopmode -halt-on-error -output-directory=$(B)

../PTAPI-v2.pdf: documento.tex guia.cpp
	mkdir -p $(B)
	$(TEX) documento.tex > $(B)/passagem1.txt || { tail -30 $(B)/passagem1.txt; exit 1; }
	$(TEX) documento.tex > $(B)/passagem2.txt || { tail -30 $(B)/passagem2.txt; exit 1; }
	$(TEX) documento.tex > $(B)/passagem3.txt || { tail -30 $(B)/passagem3.txt; exit 1; }
	cp $(B)/documento.pdf $@

clean:
	rm -rf $(B)

.PHONY: clean
