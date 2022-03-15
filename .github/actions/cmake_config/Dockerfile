FROM nvidia/cuda:10.2-devel-ubuntu18.04

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        g++ \
        wget \
        libidn11 \
        ca-certificates \
        make \
        python3-dev \
        git \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

COPY entrypoint.sh /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
