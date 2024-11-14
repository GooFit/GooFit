FROM nvidia/cuda:12.2.2-devel-ubuntu22.04

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        g++ \
        wget \
        libidn \
        ca-certificates \
        make \
        python3-dev \
        git \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

COPY entrypoint.sh /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
