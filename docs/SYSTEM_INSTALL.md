# Installing GooFit on different systems

The docker command for each system is listed at the beginning of each command, to show setup from scratch. Ignore that line if you are not using docker.

## CentOS 7

docker run -it centos
yum install epel-release -y
yum install cmake3 git -y
yum group install -y "Development Tools"
git clone --recursive https://github.com/GooFit/GooFit.git
cd GooFit
mkdir build
cd build
cmake3 ..
make


## Alpine Linux 3.5

A truly minimal system, Alpine gives you a working Docker system under 5 MB. To get a minimal install of GooFit:

```bash
docker run -it alpine sh
apk add --no-cache make cmake g++ git libexecinfo-dev
git clone --recursive https://github.com/GooFit/GooFit.git
cd GooFit
mkdir build
cd build
cmake ..
make
```

## Ubuntu 17.04

docker run -it ubuntu
apt-get update && apt-get install 
git clone --recursive https://github.com/GooFit/GooFit.git
cd GooFit
mkdir build
cd build
cmake ..
make

