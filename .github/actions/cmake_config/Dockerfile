FROM nvidia/cuda:12.9.1-devel-ubuntu24.04

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        g++ \
        wget \
        ca-certificates \
        make \
        python3-dev \
        git \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

COPY entrypoint.sh /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
