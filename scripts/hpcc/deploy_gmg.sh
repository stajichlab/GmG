#!/bin/bash
# On the HPCC, after `pixi run build` in the checkout (~/projects/GmG): install the binary as a versioned copy and swap the
# `bin/gmg` symlink. Never overwrite a binary a running job may be using (that crashed jobs with SIGSEGV/SIGBUS on 2026-10-28).
set -eu
GR=/bigdata/stajichlab/shared/projects/GmG_runs
REPO=${GMG_REPO:-$HOME/projects/GmG}
v=$(cd $REPO && git rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)
cp $REPO/target/release/gmg $GR/bin/gmg-$v
# --version must carry the commit being deployed (a cached build once reported an older hash)
head=$(cd $REPO && git rev-parse --short=9 HEAD)
$GR/bin/gmg-$v --version | grep -q "$head" || { echo "ERROR: gmg-$v --version does not show $head; rebuild with GMG_GIT_HASH=$head pixi run build" >&2; rm -f $GR/bin/gmg-$v; exit 1; }
ln -sfn gmg-$v $GR/bin/gmg.new && mv -Tf $GR/bin/gmg.new $GR/bin/gmg
echo "deployed gmg-$v"
