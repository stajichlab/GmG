#!/bin/bash
# On the HPCC, after `pixi run build` in the checkout (~/projects/GmG): install the binary as a versioned copy and swap the
# `bin/gmg` symlink. Never overwrite a binary a running job may be using (that crashed jobs with SIGSEGV/SIGBUS on 2026-10-28).
set -eu
GR=/bigdata/stajichlab/shared/projects/GmG_runs
REPO=${GMG_REPO:-$HOME/projects/GmG}
v=$(cd $REPO && git rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)
cp $REPO/target/release/gmg $GR/bin/gmg-$v
ln -sfn gmg-$v $GR/bin/gmg.new && mv -Tf $GR/bin/gmg.new $GR/bin/gmg
echo "deployed gmg-$v"
