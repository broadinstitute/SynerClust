FROM ubuntu:latest
RUN apt-get update
RUN apt-get install -y python
RUN apt-get install -y python-pip ncbi-blast+ git
RUN pip install numpy
RUN pip install networkx==1.11
RUN git clone https://github.com/SynerClust/SynerClust.git
RUN cd SynerClust && python INSTALL.py

ENV SYNERCLUST_VERSION="1.0.0"
ENV SYNERCLUST_CO="d2deb6f0a53e2ace45550b0a19e2021cedaa53b5"

ENV SYNERCLUST_HOME /SynerClust/bin
ENV PATH=${SYNERCLUST_HOME}:${PATH}

