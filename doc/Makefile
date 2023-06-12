# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = TurboRVB manual
SOURCEDIR     = source
BUILDDIR      = build
    
# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).

# Added by Nakano:
# Note: the following commands work only when path is set.

#%: Makefile
#	$(eval TurboRVB_ROOT := $(shell which turbo-genius.sh | sed -e 's@/bin/turbo-genius.sh@@g'))
#	$(eval TurboRVB_manual_ROOT := $(shell pwd))
#	echo $(TurboRVB_ROOT)
#	echo $(TurboRVB_manual_ROOT)
#	sphinx-apidoc -f -o $(TurboRVB_manual_ROOT)/source/_sources/02developers_manual/03turbo_genius_developer/98autodoc $(TurboRVB_ROOT)/turbo_genius_dist/turbo_genius/
#	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
    
#	if ! [ -e ./source/_sources/98readmes ]; then \
#	ln -s `pwd`/readme ./source/_sources/98readmes; \
#	fi

%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
