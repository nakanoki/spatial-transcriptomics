# 解析は Python のみで完結するため、R を含まない slim イメージを使う
FROM python:3.11-slim-bookworm

ENV DEBIAN_FRONTEND=noninteractive

# uv は公式イメージから取得する（バージョンを固定するため）
COPY --from=ghcr.io/astral-sh/uv:0.12.17 /uv /uvx /bin/

# uv sync が作る venv の場所を固定して PATH に通す
ENV UV_PROJECT_ENVIRONMENT=/opt/venv \
    UV_PYTHON=python3 \
    UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy \
    PATH="/opt/venv/bin:$PATH"

WORKDIR /work

# 依存関係の定義だけ先にコピーして sync する
# （ソースだけ変更した再ビルドでは、このレイヤーがキャッシュされ依存の再取得が走らない）
COPY pyproject.toml uv.lock ./
RUN uv sync --frozen --extra notebook

# ソース一式をコピー（data/ results/ は .dockerignore で除外済み）
COPY . .

# Jupyter から選択できるように ipykernel を登録する。
# venv 内（sys.prefix）に置くので、root 以外で実行しても見つかる
RUN python -m ipykernel install --sys-prefix \
      --name spatial-transcriptomics \
      --display-name "Python (spatial-transcriptomics)"

EXPOSE 8888

CMD ["jupyter", "lab", "--ip=0.0.0.0", "--no-browser", "--allow-root"]
