# syntax=docker/dockerfile:1
# gmg container image (OCI). Built with the same locked pixi environment as local development, then only the
# binary and the calibrated threshold files are copied into a slim runtime image.
#   docker build --platform linux/amd64 -t gmg .
#   docker run --rm -v "$PWD:/data" gmg guess --index /data/bfd.gmgi /data/genome.fa
# On the HPCC: apptainer pull gmg.sif docker://ghcr.io/stajichlab/gmg:<tag>
FROM ghcr.io/prefix-dev/pixi:0.83.0 AS build
WORKDIR /src
COPY pixi.toml pixi.lock ./
# solve/download the environment first so source edits do not invalidate this layer
RUN pixi install --locked
COPY . .
# the commit shown by `gmg --version` (no .git in the image build context when built from a tarball)
ARG GMG_GIT_HASH=unknown
ENV GMG_GIT_HASH=${GMG_GIT_HASH}
RUN pixi run build && ldd target/release/gmg

FROM debian:12-slim
LABEL org.opencontainers.image.title="gmg" \
      org.opencontainers.image.description="Guess My Genome: place a fungal genome or reads in the taxonomy from FracMinHash sketches" \
      org.opencontainers.image.source="https://github.com/stajichlab/GmG" \
      org.opencontainers.image.licenses="MIT"
COPY --from=build /src/target/release/gmg /usr/local/bin/gmg
COPY --from=build /src/data/thresholds /usr/share/gmg/thresholds
ENTRYPOINT ["gmg"]
CMD ["--help"]
