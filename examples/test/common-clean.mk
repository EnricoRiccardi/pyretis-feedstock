# Shared cleanup target for example tests.
#
# Each local Makefile only needs to populate:
# - CLEAN_FIND_FILES: file name patterns deleted with `find . -name ... -delete`
# - CLEAN_FIND_DIRS: directory names deleted recursively wherever they appear
# - CLEAN_RM_PATHS: root-relative paths or glob patterns removed with `rm -rf`
# - CLEAN_EXTRA_CMDS: optional plain shell commands for special cases

CLEAN_FIND_FILES ?=
CLEAN_FIND_DIRS ?=
CLEAN_RM_PATHS ?=

COMMON_CLEAN_SH := $(dir $(lastword $(MAKEFILE_LIST)))common-clean.sh

.PHONY: clean
clean:
	@CLEAN_FIND_FILES='*.pyc *.pyo $(strip $(CLEAN_FIND_FILES))' \
	CLEAN_FIND_DIRS='__pycache__ $(strip $(CLEAN_FIND_DIRS))' \
	CLEAN_RM_PATHS='$(strip $(CLEAN_RM_PATHS))' \
	CLEAN_EXTRA_CMDS='$(value CLEAN_EXTRA_CMDS)' \
	bash "$(COMMON_CLEAN_SH)"
