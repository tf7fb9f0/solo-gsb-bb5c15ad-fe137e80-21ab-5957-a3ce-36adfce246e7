FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /app

# 先单独拷贝依赖清单，利用 Docker 层缓存
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY app ./app

EXPOSE 8000

# 容器内固定监听 8000；宿主端口通过 compose 的 ports 映射（默认 8000）
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]
