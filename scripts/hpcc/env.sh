# sourced by the GmG HPCC job scripts
export GR=/bigdata/stajichlab/shared/projects/GmG_runs
export GMG=$GR/bin/gmg
export BFD_GENOMES=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes
# the git checkout on the HPCC (clone with `gh repo clone stajichlab/GmG ~/projects/GmG`; update with `git pull`)
export GMG_REPO=${GMG_REPO:-$HOME/projects/GmG}
