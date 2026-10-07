#!/usr/bin/env bash
# Prenositeľný alias pre lokálne spustenie vo WSL/Linux.
set -euo pipefail
cd "$(dirname "$0")"
exec bash start_app.sh
