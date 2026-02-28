FROM rocker/r-ver:4.4.2

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    python3 python3-pip python3-venv python3-dev \
    build-essential git curl ca-certificates \
    libssl-dev libcurl4-openssl-dev libxml2-dev \
  && rm -rf /var/lib/apt/lists/*

# Python (Jupyter + basic scientific stack)
RUN python3 -m pip install --no-cache-dir -U pip \
  && python3 -m pip install --no-cache-dir \
    jupyterlab \
    ipykernel \
    scanpy \
    squidpy \
    pyyaml \
    numpy \
    pandas \
    matplotlib

# R (IRkernel + common packages)
RUN R -q -e "install.packages(c('IRkernel','tidyverse','data.table'), repos='https://cloud.r-project.org')" \
  && R -q -e "IRkernel::installspec(user = FALSE)"

WORKDIR /work
EXPOSE 8888

CMD ["bash", "-lc", "jupyter lab --ip=0.0.0.0 --no-browser --allow-root"]
