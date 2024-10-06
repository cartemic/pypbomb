#/usr/bin/env bash
ruff check --select I --fix . && ruff check --select F401 --fix && ruff format . && ruff check