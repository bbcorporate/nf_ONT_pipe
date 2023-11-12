# starting with ngmlr
# shorthand for docker image pull
FROM ubuntu:22.04 

# first update, then install pretty generic tools
# python
# aws cli (v2)
# ncbi-blast+
# prereqs for samtools
RUN apt-get update && apt -y install make cmake ncbi-blast+ curl bzip2 unzip gawk gcc g++ python3 tar wget libncurses-dev libbz2-dev liblzma-dev zlib1g-dev git

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip awscliv2.zip
RUN ./aws/install

# use this to authenticate to gitlab? also github?
# something like this
# RUN --mount=type=ssh
# copy bifxdev repo, parts of
#COPY ../bifxdev/ONT_pipe/*.py /usr/local/bin

# bioinf stuff
##############
# ngmlr
# this release 0.2.7 has a signed integer overflow error
# **do not use**
#RUN wget https://github.com/philres/ngmlr/releases/download/v0.2.7/ngmlr-0.2.7-linux-x86_64.tar.gz && tar xzf ngmlr-0.2.7-linux-x86_64.tar.gz && mv ngmlr-0.2.7/ngmlr /usr/local/bin/
# so clone from git and compile from source
# this commit: a2a31fb
# has this fix in src/
# Fix for "Signed Integer Overflow in MAPQ"
RUN git clone https://github.com/philres/ngmlr.git && cd ngmlr && mkdir build && cd build && cmake .. && make && cd ../bin && ln -s `pwd`/ngmlr-0.2.8/ngmlr  /usr/local/bin/ngmlr
# samtools, installs into /usr/local/bin/
RUN wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 && tar -xjf samtools-1.18.tar.bz2  && cd samtools-1.18 && ./configure && make && make install 

#CMD ["ngmlr"]
#ENTRYPOINT ["/bin/bash" , "-c"]  # breaks the image when run in NF
#GATK
# has its own nf module
#snpeff
#has its own nf module, not using, hard to run
