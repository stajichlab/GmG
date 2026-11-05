#!/bin/bash
cd /bigdata/stajichlab/jstajich/projects/GmG_popgen/popgen_cocci
G=/bigdata/stajichlab/shared/projects/GmG_runs
$G/bin/gmg guess --index $G/index/bfd_s1000.gmgi sk1000/NM_7898.gsk sk1000/CA10.gsk sk1000/CA12.gsk sk1000/NM_0071.gsk 2>&1 | cut -c1-400
