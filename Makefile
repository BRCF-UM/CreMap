# CreMAP — local test environment
.PHONY: install run

install:
	Rscript scripts/install_deps.R

run:
	Rscript scripts/run_app.R
