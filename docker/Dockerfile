FROM ubuntu:latest
RUN apt-get update
RUN apt-get install -y python
RUN apt-get install -y python-pip ncbi-blast+ git
RUN pip install numpy
RUN pip install networkx==1.11
RUN git clone https://github.com/SynerClust/SynerClust.git
RUN cd SynerClust && python INSTALL.py

ENV SYNERCLUST_HOME /SynerClust/bin
ENV PATH=${SYNERCLUST_HOME}:${PATH}

