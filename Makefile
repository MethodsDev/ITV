# Makefile for managing the conda environment and package installation

VERSION=$(shell grep -oP '^__version__\s*=\s*"\K[^\"]+' src/integrative_transcriptomics_viewer/__init__.py)
VERSIONED_ENV_NAME=itv-$(VERSION)

.PHONY: install clean update

install: check-versioned-env
	mamba run -n $(VERSIONED_ENV_NAME) --no-capture-output pip install .;

check-versioned-env:
	@if mamba env list | grep -q "^$(VERSIONED_ENV_NAME)"; then \
		echo "Versioned environment '$(VERSIONED_ENV_NAME)' already exists. Skipping creation."; \
	else \
		$(MAKE) create-versioned-env; \
	fi

create-versioned-env:
	echo "Creating versioned environment '$(VERSIONED_ENV_NAME)'."
	python generate_env_yml.py
	mamba env create -f environment.yml

clean:
	@if mamba env list | grep -q "^$(VERSIONED_ENV_NAME)"; then \
		echo "Removing versioned environment '$(VERSIONED_ENV_NAME)'."; \
		mamba env remove -n $(VERSIONED_ENV_NAME); \
	else \
		echo "Versioned environment '$(VERSIONED_ENV_NAME)' does not exist. Nothing to clean."; \
	fi

update:
	@if mamba env list | grep -q "^$(VERSIONED_ENV_NAME)"; then \
		echo "Updating versioned environment '$(VERSIONED_ENV_NAME)'."; \
		mamba run -n $(VERSIONED_ENV_NAME) --no-capture-output pip install .; \
	else \
		echo "Versioned environment '$(VERSIONED_ENV_NAME)' does not exist. Run 'make install' first."; \
	fi


devinstall: check-versioned-env
	mamba run -n $(VERSIONED_ENV_NAME) --no-capture-output pip install -e .;
