#!/usr/bin/env bash
R --no-save < <(echo "devtools::document(roclets = c('rd', 'collate', 'namespace'))")
