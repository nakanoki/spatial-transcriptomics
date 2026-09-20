FROM rocker/r-ver:4.4.2

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    python3 python3-dev \
    build-essential git curl ca-certificates \
    libssl-dev libcurl4-openssl-dev libxml2-dev \
  && rm -rf /var/lib/apt/lists/*

# -----------------------
# uv
# -----------------------
COPY --from=ghcr.io/astral-sh/uv:0.12.17 /uv /uvx /bin/

# uv sync が作る venv の場所を固定して PATH に通す（旧 Dockerfile の /opt/venv を踏襲）
ENV UV_PROJECT_ENVIRONMENT=/opt/venv \
    UV_PYTHON=python3 \
    PATH="/opt/venv/bin:$PATH"

WORKDIR /work

# 依存関係の定義だけ先にコピーして sync する
# （ソースだけ変更した再ビルドでは、このレイヤーがキャッシュされ pip/uv の再取得が走らない）
COPY pyproject.toml uv.lock ./
RUN uv sync --frozen --extra notebook

# ソース一式をコピー（data/ results/ 等は .dockerignore で除外済み）
COPY . .

# ipykernel登録（JupyterでPython選択可能に）
RUN python -m ipykernel install --name python-env --display-name "Python (venv)" --user

# -----------------------
# R環境
# -----------------------
RUN R -q -e "install.packages(c('IRkernel','tidyverse','data.table'), repos='https://cloud.r-project.org')" \
  && R -q -e "IRkernel::installspec(user = FALSE)"

EXPOSE 8888

CMD ["bash", "-lc", "jupyter lab --ip=0.0.0.0 --no-browser --allow-root"]
