#!/usr/bin/env bash
set -euo pipefail

# Remove lines that exactly equal "# -*- coding: utf-8 -*-" when
# they are the first line of a Python file. If the string appears
# elsewhere in the file, print the file path (do not modify).

ROOT=${1:-.}
COOKIE="# -*- coding: utf-8 -*-"

removed_count=0
not_first_count=0

while IFS= read -r -d '' file; do
    # skip unreadable files
    [ -r "$file" ] || continue

    # read first line only and normalize CRLF
    first=$(awk 'NR==1{print; exit}' "$file" | tr -d '\r')

    if [ "$first" = "$COOKIE" ]; then
        # remove the first line safely
        tmp=$(mktemp)
        tail -n +2 "$file" > "$tmp" || true
        mv "$tmp" "$file"
        chmod --reference="$file" "$file" 2>/dev/null || true
        removed_count=$((removed_count+1))
        printf 'REMOVED %s\n' "$file"
    else
        if grep -qF "$COOKIE" "$file"; then
            printf 'NOT-FIRST %s\n' "$file"
            not_first_count=$((not_first_count+1))
        fi
    fi
done < <(find "$ROOT" -type f -name '*.py' -print0)

printf 'Done: removed=%d not-first=%d\n' "$removed_count" "$not_first_count"
