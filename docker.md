
# Docker howto

(Assuming docker is installed correctly)

```
mkdir flashpca-build
```

Put this in `flashpca-build/Dockerfile`
```
FROM ubuntu:bionic
RUN apt-get update && \
   apt-get -y install python2.7 python-pip libboost1.62-all-dev \
   libeigen3-dev git gnupg2 sudo wget ca-certificates
RUN echo 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' \
   > /etc/apt/sources.list.d/cran.list
RUN apt-key adv --keyserver keyserver.ubuntu.com \
   --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get -y install tzdata
RUN ln -fs /usr/share/zoneinfo/Australia/Melbourne /etc/localtime && \
   dpkg-reconfigure --frontend noninteractive tzdata
RUN apt-get install -y r-base r-base-dev
RUN apt-get install -y vim
RUN useradd -m flashpca-user
RUN chpasswd flashpca-user:password
WORKDIR /home/flashpca-user
USER flashpca-user
RUN wget https://github.com/yixuan/spectra/archive/v0.8.1.tar.gz && \
   tar xvf v0.8.1.tar.gz
ADD https://api.github.com/repos/gabraham/flashpca/git/refs/heads/master \
   version.json
RUN git clone https://github.com/gabraham/flashpca.git
RUN cd flashpca && \
   make all \
   EIGEN_INC=/usr/include/eigen3 \
   BOOST_INC=/usr/include/boost \
   SPECTRA_INC=$HOME/spectra-0.8.1/include &&\
   make flashpca_x86-64 \
   EIGEN_INC=/usr/include/eigen3 \
   BOOST_INC=/usr/include/boost \
   SPECTRA_INC=$HOME/spectra-0.8.1/include
```

Build the container
```
docker build flashpca-build -t 'flashpca'
```

Run the container and bind a local directory that contains your data, e.g.,
`~/Data/HapMap3`:
```
docker run -it -v ~/HapMap3:/home/flashpca-user/data flashpca
```

Run flashpca in the container:
```
flashpca-user@6dd9f366b298:~$ cd data
flashpca-user@6dd9f366b298:~/data$ ~/flashpca/flashpca --bfile chr22
[Mon May 11 11:13:30 2020] arguments: flashpca /home/flashpca-user/flashpca/flashpca --bfile chr22
[Mon May 11 11:13:30 2020] Start flashpca (version 2.1)
[Mon May 11 11:13:30 2020] blocksize: 302 (2312112 bytes per block)
[Mon May 11 11:13:30 2020] PCA begin
[Mon May 11 11:13:30 2020] PCA done
[Mon May 11 11:13:30 2020] Writing 10 eigenvalues to file eigenvalues.txt
[Mon May 11 11:13:30 2020] Writing 10 eigenvectors to file eigenvectors.txt
[Mon May 11 11:13:30 2020] Writing 10 PCs to file pcs.txt
[Mon May 11 11:13:31 2020] Writing 10 proportion variance explained to file pve.txt
[Mon May 11 11:13:31 2020] Goodbye!
```

