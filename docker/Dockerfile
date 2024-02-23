FROM python:3.9-slim

RUN apt-get update && apt-get install -y \
    build-essential \
    bzip2 \
    ca-certificates \
    apt-utils \
    zlib1g-dev \
    libssl-dev \
    libncurses5-dev \
    libsqlite3-dev \
    libreadline-dev \
    libtk8.6 \
    libgdm-dev \
    libdb4o-cil-dev \
    libpcap-dev \
    wget \
    git

COPY requirements.txt /home/
RUN python3.9 -m pip install --upgrade pip
RUN python3.9 -m pip install -r /home/requirements.txt

RUN cd /home/ && \
    git clone https://github.com/FinucaneLab/fine-mapping-inf.git || (cd fine-mapping-inf && git pull) && \
    cd /home/fine-mapping-inf/susieinf && \
    python3.9 setup.py bdist_wheel && \
    python3.9 -m pip install .

RUN cd /home/fine-mapping-inf/finemapinf && \
    python3.9 setup.py bdist_wheel && \
    python3.9 -m pip install .

VOLUME ["/root/.config"]
CMD [ "/bin/bash" ]
