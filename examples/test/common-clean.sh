#!/usr/bin/env bash
set -euo pipefail

find_files=()
find_dirs=()
rm_paths=()

if [[ -n "${CLEAN_FIND_FILES:-}" ]]; then
    read -r -a find_files <<< "${CLEAN_FIND_FILES}"
fi

if [[ -n "${CLEAN_FIND_DIRS:-}" ]]; then
    read -r -a find_dirs <<< "${CLEAN_FIND_DIRS}"
fi

if [[ -n "${CLEAN_RM_PATHS:-}" ]]; then
    read -r -a rm_paths <<< "${CLEAN_RM_PATHS}"
fi

for pattern in "${find_files[@]}"; do
    find . -name "${pattern}" -delete
done

for dirname in "${find_dirs[@]}"; do
    find . -name "${dirname}" -type d -exec rm -rf {} + 2>/dev/null || true
done

for path in "${rm_paths[@]}"; do
    rm -rf -- ${path}
done

if [[ -n "${CLEAN_EXTRA_CMDS:-}" ]]; then
    eval "${CLEAN_EXTRA_CMDS}"
fi
