set -eu

LIB=boost_1_76_0
cd pubtcrs/bin
wget https://boostorg.jfrog.io/artifactory/main/release/1.76.0/source/${LIB}.tar.gz
gunzip ${LIB}.tar.gz && tar -xf ${LIB}.tar
cd ../ && make BOOSTDIR=bin/${LIB}

