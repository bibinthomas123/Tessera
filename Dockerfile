FROM python:3.11-slim AS build
WORKDIR /src
COPY pyproject.toml README.md ./
COPY src ./src
RUN pip install --no-cache-dir build && python -m build --wheel -o /dist

FROM python:3.11-slim
RUN useradd --create-home --uid 10001 tessera
COPY --from=build /dist/*.whl /tmp/
RUN pip install --no-cache-dir "$(ls /tmp/*.whl)[service,postgres]" && rm /tmp/*.whl
USER tessera
WORKDIR /home/tessera
ENV TESSERA_DATA_DIR=/data TESSERA_CATALOG_DIRS=/catalogs TESSERA_IERS_AUTO_DOWNLOAD=0
EXPOSE 8000
CMD ["tessera", "serve", "--host", "0.0.0.0", "--port", "8000"]
