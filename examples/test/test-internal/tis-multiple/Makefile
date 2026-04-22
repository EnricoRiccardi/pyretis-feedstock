CLEAN_FIND_FILES := *out.rst* *.restart *.log*
CLEAN_RM_PATHS := tis-0* 0* report EXIT

define CLEAN_EXTRA_CMDS
	find results -name "*.html" ! -name "tis-multiple_report.html" -delete
endef

include $(shell git rev-parse --show-toplevel)/examples/test/common-clean.mk
