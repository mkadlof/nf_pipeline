FROM debian:12.4-slim

RUN apt update && \
    apt install -y --no-install-recommends openjdk-17-jre-headless \
                                           git \
                                           python3 \
                                           python3-biopython \
                                           python3-magic  \
                                           python3-pysam \
                                           curl \
                                           unzip \
                                           samtools \
                                           bwa \
                                           procps \
                                           tabix \
                                           bcftools

RUN ln -s /usr/bin/python3 /usr/bin/python

WORKDIR /opt/picard/
RUN curl -fsSL https://github.com/broadinstitute/picard/releases/download/3.1.1/picard.jar -o picard.jar

WORKDIR /opt

RUN curl -fsSL https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip -o gatk-4.5.0.0.zip && \
    unzip gatk-4.5.0.0.zip && \
    rm gatk-4.5.0.0.zip && \
    ln -s /opt/gatk-4.5.0.0 /opt/gatk

RUN git clone https://github.com/mkadlof/nf_pipeline /home/nf_pipeline
WORKDIR /home/nf_pipeline
