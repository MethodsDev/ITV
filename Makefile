# Makefile for managing the conda environment and package installation

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	VERSION=$(shell /usr/bin/grep -E '^__version__\s*=\s*' src/integrative_transcriptomics_viewer/__init__.py | grep -oE '"(.+)"' | sed 's/"//g')
endif
ifeq ($(UNAME_S),Linux)
	VERSION=$(shell grep -oP '^__version__\s*=\s*"\K[^\"]+' src/integrative_transcriptomics_viewer/__init__.py)
endif
VERSIONED_ENV_NAME=itv-$(VERSION)

.PHONY: install clean update

install: clean check-versioned-env
	conda run -n $(VERSIONED_ENV_NAME) --no-capture-output pip install .;

check-versioned-env:
	@if conda env list | grep -q "^$(VERSIONED_ENV_NAME)"; then \
		echo "Versioned environment '$(VERSIONED_ENV_NAME)' already exists. Skipping creation."; \
	else \
		$(MAKE) create-versioned-env; \
	fi

create-versioned-env:
	echo "Creating versioned environment '$(VERSIONED_ENV_NAME)'."
	python generate_env_yml.py
	conda env create -f environment.yml

clean:
	@if conda env list | grep -q "^$(VERSIONED_ENV_NAME)"; then \
		echo "Removing versioned environment '$(VERSIONED_ENV_NAME)'."; \
		conda env remove -n $(VERSIONED_ENV_NAME); \
	else \
		echo "Versioned environment '$(VERSIONED_ENV_NAME)' does not exist. Nothing to clean."; \
	fi

update:
	@if conda env list | grep -q "^$(VERSIONED_ENV_NAME)"; then \
		echo "Updating versioned environment '$(VERSIONED_ENV_NAME)'."; \
		conda run -n $(VERSIONED_ENV_NAME) --no-capture-output pip install .; \
	else \
		echo "Versioned environment '$(VERSIONED_ENV_NAME)' does not exist. Run 'make install' first."; \
	fi


devinstall: clean check-versioned-env
	conda run -n $(VERSIONED_ENV_NAME) --no-capture-output pip install -e .;
