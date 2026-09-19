FROM rocker/r-ver:4.4.2

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    python3 python3-venv python3-dev \
    build-essential git curl ca-certificates \
    libssl-dev libcurl4-openssl-dev libxml2-dev \
  && rm -rf /var/lib/apt/lists/*

# -----------------------
# Python 仮想環境
# -----------------------
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

RUN pip install --no-cache-dir --upgrade pip \
  && pip install --no-cache-dir \
      jupyterlab \
      ipykernel \
      scanpy \
      squidpy \
      pyyaml \
      numpy \
      pandas \
      matplotlib

# ipykernel登録（JupyterでPython選択可能に）
RUN python -m ipykernel install --name python-env --display-name "Python (venv)" --user

# -----------------------
# R環境
# -----------------------
RUN R -q -e "install.packages(c('IRkernel','tidyverse','data.table'), repos='https://cloud.r-project.org')" \
  && R -q -e "IRkernel::installspec(user = FALSE)"

WORKDIR /work
EXPOSE 8888

CMD ["bash", "-lc", "jupyter lab --ip=0.0.0.0 --no-browser --allow-root"]