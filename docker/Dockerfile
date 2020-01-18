FROM ubuntu:latest
WORKDIR /root

ENV rMATS=rMATS.4.0.2
RUN apt-get -y update

ENV software='libblas-dev liblapack-dev libgsl0-dev gfortran samtools wget'
ENV del=''
RUN apt-get install -y $software

RUN wget https://bootstrap.pypa.io/get-pip.py
RUN python get-pip.py
RUN pip install numpy
RUN pip install scipy

RUN wget https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz
RUN tar -zxvf 2.5.3a.tar.gz
RUN mv /root/STAR-2.5.3a/bin/Linux_x86_64_static/* /usr/local/bin/

RUN wget https://sourceforge.net/projects/rnaseq-mats/files/MATS/${rMATS}.tgz
RUN tar -zxvf ${rMATS}.tgz

RUN wget https://sourceforge.net/projects/rnaseq-mats/files/MATS/testData.tgz
RUN tar -zxvf testData.tgz

RUN rm -rf 2.5.3a.tar.gz
RUN rm -rf STAR-2.5.3a
RUN rm -rf testData*
RUN rm ${rMATS}.tgz
RUN rm get-pip.py
RUN rm -rf /var/lib/apt/lists/*

RUN mv */*UCS4/* .
RUN rm -rf ${rMATS}

ENV file=rmats.py
RUN chmod 777 ${file}
RUN ln -s ${file} RNASeq-MATS.py

RUN apt-get -y update
RUN apt-get -y install libgfortran3
RUN apt-get -y update
RUN apt-get -y upgrade
RUN apt-get -y install libgsl-dev libgsl23 libgslcblas0
RUN ln /usr/lib/x86_64-linux-gnu/libgsl.so.23 /usr/lib/x86_64-linux-gnu/libgsl.so.0

ADD ./sampleCountsSave.sh /root/

ENV PATH="${PATH}:/root"
