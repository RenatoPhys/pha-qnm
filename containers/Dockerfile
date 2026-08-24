FROM ubuntu:24.04

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    build-essential cmake ninja-build git ca-certificates \
    libcurl4-openssl-dev libssl-dev libhdf5-dev liblapack-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /work/pha-qnm
COPY . .
RUN cmake --preset release-portable && cmake --build --preset release-portable --parallel 2 \
    && ctest --preset release-portable --output-on-failure
ENTRYPOINT ["/work/pha-qnm/build/release-portable/pha-qnm"]

