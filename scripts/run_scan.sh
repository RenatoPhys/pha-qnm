#!/usr/bin/env sh
set -eu

: "${PHA_QNM_BIN:?Set PHA_QNM_BIN to the compiled executable}"
exec "$PHA_QNM_BIN" qnm scan "$@"

